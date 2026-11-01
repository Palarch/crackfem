"""Build the standalone ClaimBridge background note from its Markdown source."""
from __future__ import annotations

import re
from pathlib import Path

from docx import Document
from docx.enum.table import WD_CELL_VERTICAL_ALIGNMENT, WD_TABLE_ALIGNMENT
from docx.enum.text import WD_ALIGN_PARAGRAPH
from docx.oxml import OxmlElement
from docx.oxml.ns import qn
from docx.shared import Inches, Mm, Pt, RGBColor


ROOT = Path(__file__).resolve().parents[1]
SOURCE = ROOT / "docs" / "CLAIMBRIDGE_BACKGROUND.md"
OUTPUT = ROOT / "report" / "ClaimBridge_Background_and_FEM_Interface.docx"


def shade(cell, fill: str) -> None:
    pr = cell._tc.get_or_add_tcPr()
    el = OxmlElement("w:shd")
    el.set(qn("w:fill"), fill)
    pr.append(el)


def margins(cell, value=100) -> None:
    pr = cell._tc.get_or_add_tcPr()
    mar = OxmlElement("w:tcMar")
    for side in ("top", "start", "bottom", "end"):
        node = OxmlElement(f"w:{side}")
        node.set(qn("w:w"), str(value))
        node.set(qn("w:type"), "dxa")
        mar.append(node)
    pr.append(mar)


def keep(paragraph, with_next=False) -> None:
    pr = paragraph._p.get_or_add_pPr()
    pr.append(OxmlElement("w:keepLines"))
    if with_next:
        pr.append(OxmlElement("w:keepNext"))


def add_inline(paragraph, text: str) -> None:
    """Small Markdown subset: bold, italic, code and bare DOI URLs."""
    token = re.compile(r"(\*\*.+?\*\*|`.+?`|\*.+?\*|https://\S+)")
    pos = 0
    for match in token.finditer(text):
        if match.start() > pos:
            paragraph.add_run(text[pos:match.start()])
        value = match.group(0)
        if value.startswith("**"):
            run = paragraph.add_run(value[2:-2])
            run.bold = True
        elif value.startswith("`"):
            run = paragraph.add_run(value[1:-1])
            run.font.name = "Consolas"
            run.font.size = Pt(8.5)
        elif value.startswith("*"):
            run = paragraph.add_run(value[1:-1])
            run.italic = True
        else:
            rid = paragraph.part.relate_to(
                value.rstrip("."),
                "http://schemas.openxmlformats.org/officeDocument/2006/relationships/hyperlink",
                is_external=True,
            )
            link = OxmlElement("w:hyperlink")
            link.set(qn("r:id"), rid)
            run = OxmlElement("w:r")
            props = OxmlElement("w:rPr")
            color = OxmlElement("w:color")
            color.set(qn("w:val"), "1F4E79")
            props.append(color)
            underline = OxmlElement("w:u")
            underline.set(qn("w:val"), "single")
            props.append(underline)
            run.append(props)
            node = OxmlElement("w:t")
            node.text = value.rstrip(".")
            run.append(node)
            link.append(run)
            paragraph._p.append(link)
            if value.endswith("."):
                paragraph.add_run(".")
        pos = match.end()
    if pos < len(text):
        paragraph.add_run(text[pos:])


def add_table(doc: Document, lines: list[str]) -> None:
    parsed = [[c.strip() for c in line.strip().strip("|").split("|")] for line in lines]
    headers, rows = parsed[0], parsed[2:]
    table = doc.add_table(rows=1, cols=len(headers))
    table.alignment = WD_TABLE_ALIGNMENT.CENTER
    table.autofit = True
    for i, value in enumerate(headers):
        cell = table.rows[0].cells[i]
        cell.text = value
        shade(cell, "1F4E79")
        margins(cell)
        cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
        for p in cell.paragraphs:
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            for run in p.runs:
                run.bold = True
                run.font.color.rgb = RGBColor(255, 255, 255)
                run.font.size = Pt(8.2)
    header_pr = table.rows[0]._tr.get_or_add_trPr()
    repeat = OxmlElement("w:tblHeader")
    repeat.set(qn("w:val"), "true")
    header_pr.append(repeat)
    for row_index, values in enumerate(rows):
        row = table.add_row()
        for i, value in enumerate(values):
            cell = row.cells[i]
            cell.text = value
            margins(cell)
            cell.vertical_alignment = WD_CELL_VERTICAL_ALIGNMENT.CENTER
            if row_index % 2:
                shade(cell, "EAF2F8")
            for p in cell.paragraphs:
                p.alignment = WD_ALIGN_PARAGRAPH.LEFT if i == 0 or len(value) > 18 else WD_ALIGN_PARAGRAPH.CENTER
                p.paragraph_format.space_after = Pt(0)
                for run in p.runs:
                    run.font.size = Pt(8.1)
    doc.add_paragraph().paragraph_format.space_after = Pt(0)


def build() -> None:
    lines = SOURCE.read_text(encoding="utf-8").splitlines()
    doc = Document()
    section = doc.sections[0]
    section.page_width = Mm(210)
    section.page_height = Mm(297)
    section.top_margin = Inches(0.68)
    section.bottom_margin = Inches(0.65)
    section.left_margin = Inches(0.78)
    section.right_margin = Inches(0.78)

    normal = doc.styles["Normal"]
    normal.font.name = "Aptos"
    normal.font.size = Pt(10)
    normal._element.rPr.rFonts.set(qn("w:eastAsia"), "Aptos")
    normal.paragraph_format.space_after = Pt(6)
    normal.paragraph_format.line_spacing = 1.08
    title_style = doc.styles["Title"]
    title_style.font.name = "Aptos Display"
    title_style.font.size = Pt(24)
    title_style.font.bold = True
    title_style.font.color.rgb = RGBColor(0, 0, 0)
    for name, size in (("Heading 1", 15), ("Heading 2", 12)):
        style = doc.styles[name]
        style.font.name = "Aptos Display"
        style.font.size = Pt(size)
        style.font.bold = True
        style.font.color.rgb = RGBColor(0, 0, 0)
        style.paragraph_format.space_before = Pt(10)
        style.paragraph_format.space_after = Pt(4)
    caption_style = doc.styles["Caption"]
    caption_style.font.name = "Aptos"
    caption_style.font.size = Pt(8.5)
    caption_style.font.italic = True
    caption_style.font.color.rgb = RGBColor(70, 70, 70)

    title = doc.add_paragraph(style="Title")
    title.alignment = WD_ALIGN_PARAGRAPH.CENTER
    title.add_run(lines[0][2:])
    for element in (title._p.get_or_add_pPr(), title_style._element.get_or_add_pPr()):
        border = element.find(qn("w:pBdr"))
        if border is not None:
            element.remove(border)
    sub = doc.add_paragraph()
    sub.alignment = WD_ALIGN_PARAGRAPH.CENTER
    run = sub.add_run("Supplement to the FEM course project")
    run.bold = True
    run.font.size = Pt(12)
    date = doc.add_paragraph("Prepared for supervisor review · 18 October 2026")
    date.alignment = WD_ALIGN_PARAGRAPH.CENTER
    date.runs[0].font.size = Pt(9)
    date.runs[0].font.color.rgb = RGBColor(90, 90, 90)

    i = 1
    in_references = False
    while i < len(lines):
        line = lines[i].strip()
        if not line:
            i += 1
            continue
        if line.startswith("## "):
            p = doc.add_heading(line[3:], level=1)
            keep(p, True)
            if line == "## 3 Interface contract":
                p.paragraph_format.page_break_before = True
            in_references = line == "## References"
        elif line.startswith("!["):
            match = re.match(r"!\[(.*?)\]\((.*?)\)", line)
            if not match:
                raise ValueError(f"bad image line: {line}")
            image_path = (SOURCE.parent / match.group(2)).resolve()
            p = doc.add_paragraph()
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            p.paragraph_format.keep_with_next = True
            width = Inches(2.25) if "case2_" in image_path.name else Inches(6.65)
            p.add_run().add_picture(str(image_path), width=width)
        elif line.startswith("|"):
            block = []
            while i < len(lines) and lines[i].strip().startswith("|"):
                block.append(lines[i].strip())
                i += 1
            add_table(doc, block)
            continue
        elif re.match(r"^\d+\. ", line):
            if in_references:
                p = doc.add_paragraph()
                p.paragraph_format.left_indent = Inches(0.18)
                p.paragraph_format.first_line_indent = Inches(-0.18)
                p.paragraph_format.space_after = Pt(3)
                p.paragraph_format.line_spacing = 1.0
                add_inline(p, line)
                for run in p.runs:
                    run.font.size = Pt(8.5)
            else:
                p = doc.add_paragraph(style="List Number")
                add_inline(p, re.sub(r"^\d+\. ", "", line))
                p.paragraph_format.space_after = Pt(4)
        elif line.startswith("- "):
            p = doc.add_paragraph(style="List Bullet")
            add_inline(p, line[2:])
            p.paragraph_format.space_after = Pt(4)
        elif line.startswith("> "):
            p = doc.add_paragraph()
            p.paragraph_format.left_indent = Inches(0.35)
            p.paragraph_format.right_indent = Inches(0.35)
            p.paragraph_format.space_before = Pt(5)
            p.paragraph_format.space_after = Pt(8)
            r = p.add_run(line[2:])
            r.italic = True
        elif line.startswith("*") and line.endswith("*") and not line.startswith("**"):
            p = doc.add_paragraph(style="Caption")
            p.alignment = WD_ALIGN_PARAGRAPH.CENTER
            add_inline(p, line)
            keep(p)
        else:
            p = doc.add_paragraph()
            add_inline(p, line)
            keep(p, line.endswith(":"))
        i += 1

    header = section.header.paragraphs[0]
    header.text = "ClaimBridge Background and FEM Interface Note"
    header.alignment = WD_ALIGN_PARAGRAPH.RIGHT
    header.runs[0].font.size = Pt(8)
    header.runs[0].font.color.rgb = RGBColor(110, 110, 110)
    footer = section.footer.paragraphs[0]
    footer.alignment = WD_ALIGN_PARAGRAPH.CENTER
    footer.add_run("FEM course project · supervisor background note · ")
    page = OxmlElement("w:fldSimple")
    page.set(qn("w:instr"), "PAGE")
    footer._p.append(page)
    for run in footer.runs:
        run.font.size = Pt(8)
        run.font.color.rgb = RGBColor(110, 110, 110)

    doc.core_properties.title = "ClaimBridge Background and FEM Interface Note"
    doc.core_properties.subject = "Background note for the ClaimBridge measurement interface used by the FEM course project"
    doc.core_properties.comments = "Prepared from committed project evidence; field validation remains pending."
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    doc.save(OUTPUT)
    print(f"wrote {OUTPUT.relative_to(ROOT)}")


if __name__ == "__main__":
    build()

