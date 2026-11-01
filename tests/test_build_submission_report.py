"""The markdown subset that build_submission_report.py renders into the Word report."""
import importlib.util
from pathlib import Path

import pytest

docx = pytest.importorskip("docx")          # packaging dependency only, not a project dependency

_spec = importlib.util.spec_from_file_location(
    "build_submission_report", Path(__file__).resolve().parents[1] / "scripts" / "build_submission_report.py")
bsr = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(bsr)

ROOT = Path(__file__).resolve().parents[1]


def test_section_extracts_one_heading_and_its_body():
    md = ["# 7 Old", "x", "# 8 New block", "## 8.1 Part", "text", "# 9 Next", "y"]
    assert bsr.section(md, "# 8 New block") == ["# 8 New block", "## 8.1 Part", "text"]


def test_block_runs_from_its_first_line_to_the_next_heading():
    md = ["## 1.3 Read", "intro", "**Progress at a glance.** rows", "| A |", "|---|", "| 1 |", "", "## 1.4 Next"]
    assert bsr.block(md, "**Progress at a glance.**") == md[2:7]


def test_bold_leads_splits_a_plain_label_off_its_paragraph():
    doc = docx.Document()
    p = doc.add_paragraph("Found: five problems.")
    done = doc.add_paragraph()
    done.add_run("Planned:").bold = True
    done.add_run(" a design.")
    bsr.bold_leads(doc)
    assert [(r.text, r.bold) for r in p.runs] == [("Found:", True), (" five problems.", None)]
    assert [(r.text, r.bold) for r in done.runs] == [("Planned:", True), (" a design.", None)]


def test_tidy_tables_sizes_equal_columns_and_keeps_rows_whole():
    doc = docx.Document()
    table = doc.add_table(rows=2, cols=2)
    for row, texts in zip(table.rows, (("#", "Problem found"), ("P1", "a long description of the problem found"))):
        for cell, text in zip(row.cells, texts):
            cell.text = text
            cell.width = bsr.Inches(3.23)
    bsr.tidy_tables(doc)
    w = [table.rows[0].cells[c].width for c in range(2)]
    assert w[1] > 2 * w[0]
    for r, row in enumerate(table.rows):
        tr_pr = row._tr.trPr
        assert tr_pr.find(bsr.qn("w:cantSplit")) is not None
        assert (tr_pr.find(bsr.qn("w:tblHeader")) is not None) == (r == 0)


def test_render_writes_headings_lists_tables_and_figures_in_order():
    doc = docx.Document()
    anchor = doc.add_paragraph("END")
    md = ["# 8 Block", "", "**Planned:** close `the gaps`.", "", "- one", "- two", "", "1. first",
          "", "| A | B |", "|---|---|", "| 1 | 2 |", "",
          "![Figure 11. A caption.](../results/figures/fig0_workflow.png)"]
    bsr.render(doc, anchor._p, md, ROOT / "report")
    texts = [p.text for p in doc.paragraphs]
    assert texts[:6] == ["8 Block", "Planned: close the gaps.", "one", "two", "1.\tfirst", ""]
    assert doc.paragraphs[0].style.name == "Heading 1"
    assert doc.paragraphs[2].style.name == "List Bullet"
    lead = doc.paragraphs[1].runs
    assert lead[0].bold and lead[0].text == "Planned:" and lead[2].font.name == "Consolas"
    assert doc.paragraphs[6].text == "Figure 11. A caption." and doc.paragraphs[6].runs[0].italic
    assert texts[-1] == "END"
    (table,) = doc.tables
    assert [c.text for c in table.rows[0].cells] == ["A", "B"] and table.rows[0].cells[0].paragraphs[0].runs[0].bold
    body = list(doc.element.body.iterchildren())
    assert body.index(table._tbl) < body.index(doc.paragraphs[5]._p)    # table sits before the figure


def test_code_runs_match_the_draft_size_and_columns_follow_their_text():
    doc = docx.Document()
    anchor = doc.add_paragraph("END")
    md = ["Run `run_all.py` now.", "", "| A | B |", "|---|---|",
          "| x | a much longer description that needs most of the width of the page |"]
    bsr.render(doc, anchor._p, md, ROOT / "report")
    assert doc.paragraphs[0].runs[1].font.size == bsr.CODE_SIZE
    (table,) = doc.tables
    w = [table.rows[0].cells[c].width for c in range(2)]
    assert w[1] > 2 * w[0]
    assert sum(w) == pytest.approx(bsr.TEXT_WIDTH, rel=0.01)
