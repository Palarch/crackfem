"""Build the Word submission report from the reviewed draft.

This packaging script uses the existing formatted DOCX as its template. Run it with a Python
environment that provides python-docx (the numerical project itself does not depend on it).
"""
from __future__ import annotations

import copy
import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.text.paragraph import Paragraph
from docx.text.run import Run
from docx.shared import Emu, Inches, Pt


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "report" / "crackfem_report_DRAFT.docx"
OUTPUT = ROOT / "report" / "crackfem_report_FINAL.docx"
WORKFLOW = ROOT / "results" / "figures" / "fig0_workflow.png"
MARKDOWN = ROOT / "report" / "report_draft.md"
TESTS = "424"                 # python -m pytest -q, 1 November 2026, with the private pipeline
TESTS_DATE = "1 November 2026"
REPOSITORY = "https://github.com/Palarch/crackfem"
SECTION_8 = "# 8 Block 6: Closing the proposal's outcomes (24 October)"
PROGRESS = "**Progress at a glance.**"
LEADS = ("Planned:", "Found:", "Changed and why:", "Numbers:")
COVER_SPACE = Pt(150)
APPENDIX_B = "# Appendix B: Hand calculation of one beam element"
RENUMBER = ("13 Student decisions", "12 Conclusions", "11 Limitations", "10 Discussion", "9 Verification summary",
            "8 Where FEM sits")


def find_start(doc: Document, prefix: str) -> Paragraph:
    matches = [p for p in doc.paragraphs if p.text.startswith(prefix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one paragraph starting {prefix!r}, found {len(matches)}")
    return matches[0]


def replace_start(doc: Document, prefix: str, text: str) -> Paragraph:
    p = find_start(doc, prefix)
    p.text = text
    return p


def insert_after(paragraph: Paragraph, text: str = "", style: str | None = None) -> Paragraph:
    element = OxmlElement("w:p")
    paragraph._p.addnext(element)
    new = Paragraph(element, paragraph._parent)
    if style:
        new.style = style
    if text:
        new.add_run(text)
    return new


# Markdown subset of report_draft.md rendered into the draft's Word styles, so new sections are
# written once (in the markdown) and match figures 1-10 and the tables of the draft.
FIG_WIDTH = Inches(4.92)            # width of figures 1-10 in the draft
CAPTION_SIZE = Emu(120650)          # 9.5 pt italic captions
TABLE_SIZE = Pt(9)
CODE_SIZE = Emu(120650)             # 9.5 pt Consolas in the draft's body text
TEXT_WIDTH = Inches(6.46)           # A4 less the draft's 0.906 in margins
MIN_COLUMN = Inches(0.75)
NUMBERED_INDENT = (Emu(269875), Emu(-215900))   # "1.<tab>" items of Section 12
INLINE = re.compile(r"(\*\*.+?\*\*|`.+?`)")
IMAGE = re.compile(r"^!\[(.+)\]\((.+)\)$")
NUMBERED = re.compile(r"^(\d+)\. (.+)$")


def section(md: list[str], heading: str) -> list[str]:
    """The lines of one top-level section: its heading up to the next '# ' heading."""
    start = md.index(heading)
    end = next((i for i in range(start + 1, len(md)) if md[i].startswith("# ")), len(md))
    return md[start:end]


def block(md: list[str], first: str) -> list[str]:
    """The lines from the one starting with first up to the next heading of any level."""
    start = next(i for i, line in enumerate(md) if line.startswith(first))
    end = next((i for i in range(start + 1, len(md)) if md[i].startswith("#")), len(md))
    return md[start:end]


def add_runs(paragraph: Paragraph, text: str, size=None) -> None:
    for token in INLINE.split(text):
        if not token:
            continue
        if token.startswith("**"):
            run = paragraph.add_run(token[2:-2])
            run.bold = True
        elif token.startswith("`"):
            run = paragraph.add_run(token[1:-1])
            run.font.name = "Consolas"
            run.font.size = CODE_SIZE
        else:
            run = paragraph.add_run(token)
        if size is not None:
            run.font.size = size


def column_widths(rows: list[list[str]]) -> list:
    """Column widths proportional to the square root of the mean text length of each column (so a
    short status column is not squeezed), at least MIN_COLUMN, filling the text width."""
    weight = [(sum(len(r[c]) for r in rows) / len(rows)) ** 0.5 for c in range(len(rows[0]))]
    share = [max(m / sum(weight) * TEXT_WIDTH, MIN_COLUMN) for m in weight]
    return [Emu(int(w * TEXT_WIDTH / sum(share))) for w in share]


def render(doc: Document, anchor, md: list[str], base: Path) -> None:
    """Insert the markdown lines before the XML element anchor: '#'/'##' headings, paragraphs
    with **bold** and `code`, '- ' bullets, '1. ' numbered items, pipe tables and
    ![caption](figure) figures (paths relative to base)."""
    def paragraph(style=None) -> Paragraph:
        element = OxmlElement("w:p")
        anchor.addprevious(element)
        p = Paragraph(element, doc._body)
        if style:
            p.style = style
        return p

    i = 0
    while i < len(md):
        line = md[i].rstrip()
        if not line:
            i += 1
        elif line.startswith("|"):
            rows = []
            while i < len(md) and md[i].startswith("|"):
                cells = [c.strip() for c in md[i].strip().strip("|").split("|")]
                if not all(set(c) <= set("-: ") for c in cells):
                    rows.append(cells)
                i += 1
            table = doc.add_table(rows=len(rows), cols=len(rows[0]))
            table.style = "Table Grid"
            table.alignment = WD_TABLE_ALIGNMENT.CENTER
            table.autofit = False
            widths = column_widths(rows)
            for r, cells in enumerate(rows):
                for c, text in enumerate(cells):
                    table.rows[r].cells[c].width = widths[c]
                    cell_p = table.rows[r].cells[c].paragraphs[0]
                    add_runs(cell_p, f"**{text}**" if r == 0 and text else text, TABLE_SIZE)
            anchor.addprevious(table._tbl)
        elif IMAGE.match(line):
            caption, rel = IMAGE.match(line).groups()
            p = paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.keep_with_next = True
            p.add_run().add_picture(str((base / rel).resolve()), width=FIG_WIDTH)
            cap = paragraph()
            cap.alignment = WD_ALIGN_PARAGRAPH.CENTER
            run = cap.add_run(caption)
            run.italic, run.font.size = True, CAPTION_SIZE
            i += 1
        elif line.startswith("## ") or line.startswith("# "):
            level = "Heading 2" if line.startswith("## ") else "Heading 1"
            paragraph(level).add_run(line.split(" ", 1)[1])
            i += 1
        elif line.startswith("- "):
            add_runs(paragraph("List Bullet"), line[2:])
            i += 1
        elif NUMBERED.match(line):
            number, text = NUMBERED.match(line).groups()
            p = paragraph()
            p.paragraph_format.left_indent, p.paragraph_format.first_line_indent = NUMBERED_INDENT
            add_runs(p, f"{number}.	{text}")
            i += 1
        else:
            add_runs(paragraph(), line)
            i += 1


def bold_leads(doc: Document) -> None:
    """Block labels ('Planned:', 'Found:', ...) in bold, also where replace_start rewrote the
    paragraph as one plain run."""
    for p in doc.paragraphs:
        lead = next((label for label in LEADS if p.text.startswith(label)), None)
        if lead is None or p.runs[0].bold or not p.runs[0].text.startswith(lead):
            continue
        first = p.runs[0]
        if first.text != lead:
            tail = Run(copy.deepcopy(first._r), p)
            first._r.addnext(tail._r)
            tail.text = first.text[len(lead):]
        first.text, first.bold = lead, True


def _row_flag(row, tag: str) -> None:
    """Add w:cantSplit or w:tblHeader to a row in schema order (both come around w:trHeight)."""
    tr_pr = row._tr.get_or_add_trPr()
    if tr_pr.find(qn(tag)) is not None:
        return
    after = {"w:cantSplit": ("w:trHeight", "w:tblHeader", "w:tblCellSpacing", "w:jc", "w:hidden"),
             "w:tblHeader": ("w:tblCellSpacing", "w:jc", "w:hidden")}[tag]
    flag = OxmlElement(tag)
    successor = next((child for child in tr_pr if child.tag in {qn(t) for t in after}), None)
    if successor is None:
        tr_pr.append(flag)
    else:
        successor.addprevious(flag)


def tidy_tables(doc: Document) -> None:
    """Columns sized by their text where the draft left them equal, no row split across a page,
    and the header row repeated after a page break and never left alone at the foot of a page."""
    for table in doc.tables:
        if len({cell.width for cell in table.rows[0].cells}) == 1:
            widths = column_widths([[cell.text for cell in row.cells] for row in table.rows])
            table.autofit = False
            for row in table.rows:
                for cell, width in zip(row.cells, widths):
                    cell.width = width
        for r, row in enumerate(table.rows):
            _row_flag(row, "w:cantSplit")
            if r == 0:
                _row_flag(row, "w:tblHeader")
                for cell in row.cells:
                    for paragraph in cell.paragraphs:
                        paragraph.paragraph_format.keep_with_next = True


def add_cover_and_contents(doc: Document) -> None:
    """The title block alone on page 1 (no page number), then a contents page whose TOC field
    Word fills when the PDF is exported, then the Summary on a new page."""
    title = doc.paragraphs[0]
    if title.style.name != "Title":
        raise RuntimeError("the draft must open with its Title paragraph")
    title.paragraph_format.space_before = COVER_SPACE
    doc.sections[0].different_first_page_header_footer = True
    summary = find_start(doc, "Summary")
    heading = Paragraph(OxmlElement("w:p"), summary._parent)
    summary._p.addprevious(heading._p)
    heading.style = "TOC Heading"
    heading.paragraph_format.page_break_before = True
    heading.add_run("Contents")
    toc = Paragraph(OxmlElement("w:p"), summary._parent)
    summary._p.addprevious(toc._p)
    for kind, text in (("begin", None), ("instr", 'TOC \\o "1-2" \\h \\z \\u'), ("separate", None),
                       ("text", "Contents not yet generated."), ("end", None)):
        run = toc.add_run()
        if kind == "text":
            run.text = text
        elif kind == "instr":
            instr = OxmlElement("w:instrText")
            instr.set(qn("xml:space"), "preserve")
            instr.text = text
            run._r.append(instr)
        else:
            char = OxmlElement("w:fldChar")
            char.set(qn("w:fldCharType"), kind)
            run._r.append(char)
    summary.paragraph_format.page_break_before = True


def all_text(doc: Document) -> str:
    parts = [p.text for p in doc.paragraphs]
    for table in doc.tables:
        for row in table.rows:
            parts.extend(cell.text for cell in row.cells)
    for section in doc.sections:
        parts.extend(p.text for p in section.header.paragraphs)
        parts.extend(p.text for p in section.footer.paragraphs)
    return "\n".join(parts)


def add_block6(doc: Document) -> None:
    """Block 6 (24 October): new Section 8 and Appendix B rendered from the markdown, old Sections
    8-13 renumbered 9-14, and the summary, Section 2, verification table, conclusions and
    Appendix A updated to match report_draft.md."""
    md = MARKDOWN.read_text(encoding="utf-8").splitlines()
    for old in RENUMBER:
        number, rest = old.split(" ", 1)
        heading = find_start(doc, old)
        heading.text = f"{int(number) + 1} {heading.text.split(' ', 1)[1]}"
    render(doc, find_start(doc, "9 Where FEM sits")._p, section(md, SECTION_8), MARKDOWN.parent)
    render(doc, doc.element.body.sectPr, section(md, APPENDIX_B), MARKDOWN.parent)

    insert_after(
        find_start(doc, "The structural model matters more than the camera"),
        "Putting all the damage in one element, as the proposal planned, makes the answer depend on the mesh: "
        "at η = 0.36 one element gives 27.1 mm of truck deflection at N = 10 (inside L/800) and 14.2 mm at "
        "N = 40, where the model with a cracked zone set by the moments gives 72.1 mm at every mesh.",
        "List Bullet",
    )
    insert_after(
        find_start(doc, "a ±2·GSD band on the result"),
        "Section 8.4 maps each of these items and the proposal's five Expected Outcomes to the place where it "
        "was delivered.",
    )
    table = next(t for t in doc.tables if t.rows[-1].cells[0].text == "Whole code")
    last = table.rows[-1]
    for cells in (("Hand calculation, one element (Appendix B)",
                   "12EI/L³ … 4EI/L, qL/2 and qL²/12 by hand; 5qL⁴/384EI; qL/2",
                   "equal to the code; deflection and reactions to 1e-9"),
                  ("Damage option of the Stage-1 girder",
                   "damage everywhere = distributed model; no damage = linear uncracked solution",
                   "to 1e-12 and 1e-9")):
        row = table.add_row()
        for cell, text in zip(row.cells, cells):
            add_runs(cell.paragraphs[0], text, TABLE_SIZE)
        last._tr.addprevious(row._tr)
    for run in last.cells[1].paragraphs[0].runs:
        run.text = run.text.replace("401", TESTS)
    third = find_start(doc, "3.\tThe deck model changes the engineering verdict")
    fourth = insert_after(
        third,
        "4.\tPutting all the damage in one element, as the proposal planned, gives an answer set by the mesh "
        "(27.1 mm at N = 10, 14.2 mm at N = 40, against 72.1 mm with the cracked zone set by the moments): the "
        "damaged length has to come from the moments, not from the element size.",
    )
    fourth.paragraph_format.left_indent, fourth.paragraph_format.first_line_indent = NUMBERED_INDENT
    fifth = find_start(doc, "4.\tThe method is limited by what the camera can see")
    fifth.runs[0].text = fifth.runs[0].text.replace("4.\t", "5.\t", 1)
    stage2 = find_start(doc, "python scripts/run_stage2.py")
    bullet = insert_after(stage2, style="List Bullet")
    add_runs(bullet, "`python scripts/run_proposal_checks.py`: Block 6: `results/proposal_checks.json` and figures "
                     "15–17 (about 3 minutes; it reads stage1.json and stage2.json)")
    demo = insert_after(bullet, style="List Bullet")
    add_runs(demo, "`python scripts/build_demo_data.py`: the data of the browser demo, `docs/demo_data.js` (about "
                   "5 minutes); the demo itself is at https://palarch.github.io/crackfem/ and runs no FEM")
    both = find_start(doc, "Both scripts start from the committed")
    both.runs[0].text = both.runs[0].text.replace("Both scripts", "The first two scripts", 1)
    if not (fifth.text.startswith("5.\t") and both.text.startswith("The first two scripts")
            and TESTS in last.cells[1].text):
        raise RuntimeError("a Block 6 run-level replacement did not apply")


def main() -> None:
    doc = Document(SOURCE)
    md = MARKDOWN.read_text(encoding="utf-8").splitlines()
    replace_start(doc, "Repository: STUDENT: GitHub link", f"1 November 2026 · code, tests and results: {REPOSITORY}")
    replace_start(
        doc,
        "A bridge inspection rates a crack by its width",
        "ClaimBridge existed before this course project as part of the author's thesis work on image-based "
        "bridge inspection. It supplies a measured crack width and a resolution band; it does not calculate "
        "the structural consequences. The new FEM course work begins at that boundary and maps the "
        "ClaimBridge measurement to a model-conditional estimate of remaining prestress (η, a fraction of "
        "the design value), stiffness, truck deflection and load sharing.",
    )
    replace_start(
        doc,
        "The code (Python, 401 automated tests)",
        f"The code (Python, {TESTS} automated tests) was written by Claude (Anthropic's AI). The student's final "
        "code and oral review is recorded separately before submission; Section 1.4 states the arrangement.",
    )
    replace_start(
        doc,
        "The rule-engine rating is 5",
        "The rule-engine rating is 5 ‘Excellent’ for every measured case, including a crack that the stand-in "
        "model maps to a 79% prestress loss: the rating carries no information about prestress, which is the "
        "gap the indicator addresses.",
    )
    replace_start(
        doc,
        "Thai highway bridges are inspected visually",
        "Thai highway bridges are inspected visually and each defect gets a condition rating. A crack is "
        "rated by its width. In a reinforced concrete beam a hairline flexural crack is expected and may be "
        "acceptable. In a prestressed girder designed to remain in compression under permanent load, an open "
        "flexural crack at rest is evidence that the local precompression may have been exhausted. Under the "
        "adopted model, and provided the opening is not residual damage from a past overload, its width carries "
        "information about the prestress that remains; a categorical rating does not use that information.",
    )
    replace_start(doc, "1.2 The AI + FEM idea", "1.2 From ClaimBridge background to the FEM deliverables")
    replace_start(
        doc,
        "The project joins two tools.",
        "The project follows one continuous chain: ClaimBridge thesis background → crack width w ± u → "
        "inverse FEM → remaining prestress η → nonlinear girder and grillage response → engineering outputs. "
        "The hand-over variable is the measured width and its resolution band. ClaimBridge is the upstream "
        "measurement source; the structural-analysis layer is the work submitted to the FEM course.",
    )
    replace_start(
        doc,
        "AI part (measurement). ClaimBridge",
        "AI part (measurement). ClaimBridge, the author's thesis inspection pipeline (private), receives a "
        "crack mask and surface scale, measures the crack width, and passes the result to a rule engine. The "
        "upstream image-to-mask AI remains outside this course repository and is not evaluated here; the "
        "interface tested by this project begins with the mask and GSD. The pipeline reports the width as "
        "w ± u with u = 2·GSD and refuses cracks narrower than 3 pixels.",
    )
    fem_part = replace_start(
        doc,
        "FEM part (mechanics). The finite element model",
        "FEM part (mechanics). The finite element model of the girder or the whole deck maps that width to a "
        "model-conditional estimate of the remaining prestress η (an inverse problem), and η to stiffness, "
        "deflection and load sharing (a forward nonlinear analysis).",
    )
    deliverables_intro = insert_after(fem_part, "The FEM course deliverables are:")
    cursor = deliverables_intro
    for item in (
        "a verified Stage-1 Euler–Bernoulli beam model with element assembly, equivalent nodal loads and secant stiffness;",
        "an inverse solver for η with a Monte Carlo band from the ClaimBridge width measurement;",
        "a nonlinear five-girder grillage that captures stiffness-dependent load redistribution and truck deflection; and",
        "a reproducible evidence package containing tests, result files, figures, sources and decision records.",
    ):
        cursor = insert_after(cursor, item, "List Bullet")
    insert_after(
        cursor,
        "ClaimBridge is not resubmitted as the FEM contribution; it supplies the observation that makes the "
        "structural inverse problem possible.",
    )
    flow_anchor = find_start(doc, "The result is an advisory structural indicator")
    image_p = insert_after(flow_anchor)
    image_p.alignment = WD_ALIGN_PARAGRAPH.CENTER
    image_p.paragraph_format.keep_with_next = True
    image_p.add_run().add_picture(str(WORKFLOW), width=Inches(6.7))
    caption = insert_after(
        image_p,
        "Figure 0. Boundary and data flow from the crack image to the engineering outputs. The private "
        "upstream image-to-mask AI is outside the course repository.",
        "Caption",
    )
    caption.alignment = WD_ALIGN_PARAGRAPH.CENTER

    replace_start(
        doc,
        "The course asked for a report that tells the progress as it happened.",
        "The course asked for a report that tells the progress as it happened. Sections 3 to 8 follow the "
        "work blocks in order. Each block says what was planned, what was found, what changed and why, and the "
        "numbers it produced. Every number comes from the result files in the repository "
        "(results/stage1.json, results/stage2.json) or from test output, and every non-obvious decision has a "
        "dated line in docs/DECISIONS.md. Section 8.4 maps every promise of the proposal to where it was "
        "delivered. Section 9 shows where the course's FE methods sit in the code, Section 10 collects the "
        "verification, and Sections 11 to 13 discuss the results, limitations and next steps. "
        "The separate ClaimBridge background note documents the inherited measurement side for a reader who "
        "was not involved in the thesis.",
    )
    render(doc, find_start(doc, "1.4 Authorship")._p, block(md, PROGRESS), MARKDOWN.parent)

    replace_start(
        doc,
        "The code was written by Claude",
        "The code was written by Claude (Anthropic's AI); this was the student's choice at the start of the "
        "design (DECISIONS 2026-10-07), because the final exam fell in the middle of the project. The student "
        "made the engineering decisions recorded as ‘student decision’ in the log, including staged "
        "construction and the DOH truck factor. Each block in docs/PROGRESS_LOG.md ends with a reading guide, "
        "three explain-back questions and a ‘Reviewed by student’ line. The student declared the review complete "
        "on 1 November 2026 (Section 14). The report text was drafted by Claude from the progress log and the "
        "decision journal; the decisions marked as the student's are the student's.",
    )
    replace_start(
        doc,
        "Changed and why: the unknown became the remaining prestress",
        "Changed and why: the unknown became the remaining prestress η because, under the adopted flexural "
        "model, a crack visible at rest indicates that local precompression has been exhausted; η is therefore "
        "the model quantity inferred from the width. Two other approaches were rejected: a forward-only chain "
        "from width to moment to stiffness (no real nonlinearity), and continuous spans (no evidence that the "
        "spans are continuous). The work was split into Stage 1, one girder with the whole chain, complete and "
        "submittable on its own, and Stage 2, the grillage of the deck as the main deliverable. Crack width and "
        "stiffness both use Eurocode 2; ACI Branson is kept as a cross-check. The machine-learning surrogate was "
        "dropped because one solve costs milliseconds (Section 6 gives the timing).",
    )
    replace_start(
        doc,
        "On the AI side, the pipeline provides",
        "On the AI side, the pipeline provides the measurement and its uncertainty; the FE model provides the "
        "mechanics. A second AI, Claude, wrote the code, with the student's engineering decisions recorded in "
        "the decision journal and the final student review tracked separately. The independent checks in "
        "Section 10 support the numerical results; authorship alone does not.",
    )
    replace_start(
        doc,
        "STUDENT: add your own discussion point",
        "For a bridge owner, the present output should be read as a screening indicator rather than a capacity "
        "assessment. Its most actionable finding is the measurement requirement: the field GSD of 0.38 mm/px "
        "cannot resolve the crack widths predicted by this stand-in girder, whereas close-range images near "
        "0.02–0.05 mm/px can support the inverse. Before acting on the estimated η, the owner would need the "
        "as-built geometry, the history and cause of the crack, and an independent structural assessment.",
    )
    heading = replace_start(doc, "13 Student review and reflection", "13 Student decisions and review before submission")
    review = replace_start(
        doc,
        "STUDENT: this section is yours",
        "The student made three material engineering decisions recorded during the project: staged "
        "construction was retained, the HS20-44 load was multiplied by the DOH factor 1.30, and the whole-deck "
        "grillage became the main result while the single girder remained the fallback and comparison model. "
        "These choices changed the decompression threshold, service response and interpretation of the same "
        "measured crack.",
    )
    insert_after(
        review,
        "The implementation was produced with AI assistance, so submission requires a separate oral and code "
        "review by the student. The review record is report/STUDENT_REVIEW.md. It covers the beam and grid "
        "elements, assembly and partition, the secant loop, the inverse, the uncertainty statement, and the "
        "distinction between a synthetic consistency check and field validation. The student declared the "
        "review complete on 1 November 2026; the review record states which answers were written down.",
    )
    replace_start(
        doc,
        "1.\tA crack width measured at rest on a prestressed girder can be turned",
        "1.\tUnder the adopted crack and structural models, a crack width measured at rest can be mapped to "
        "estimates of remaining prestress, stiffness, truck deflection and deck load sharing, with a band from "
        "the measurement. The rule-engine rating does not provide this mechanical interpretation.",
    )
    replace_start(
        doc,
        "python -m pytest -q: 401 tests",
        f"python -m pytest -q: {TESTS} passed with the configured private pipeline (captured {TESTS_DATE}); "
        "without that configuration the 3 live-pipeline tests are skipped",
    )
    add_block6(doc)

    for paragraph in doc.paragraphs:
        if "S01–S27" in paragraph.text:
            paragraph.text = paragraph.text.replace("S01–S27", "S01–S30")
        if "S01-S27" in paragraph.text:
            paragraph.text = paragraph.text.replace("S01-S27", "S01-S30")

    tidy_tables(doc)
    bold_leads(doc)
    add_cover_and_contents(doc)
    find_start(doc, "Appendix A").paragraph_format.page_break_before = True

    # Ensure the report itself contains no visible drafting placeholder.
    forbidden = ("STUDENT:", "[[", "]]", "Repository: STUDENT", "Draft of 2026-10-15", "TESTTOTAL", "TESTCOUNT")
    text = all_text(doc)
    leftovers = [token for token in forbidden if token in text]
    if leftovers:
        raise RuntimeError(f"draft placeholders remain: {leftovers}")
    doc.core_properties.title = "From a Measured Crack Width to Remaining Prestress Stiffness and Load Sharing"
    doc.core_properties.subject = "FEM course project progress report"
    doc.core_properties.author = "Palid Archavachalee"
    doc.save(OUTPUT)
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
