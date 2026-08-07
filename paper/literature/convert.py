#!/usr/bin/env python3
"""Convert downloaded PDFs in pdfs/ to plain text in text/ using pypdf."""
import sys
from pathlib import Path
from pypdf import PdfReader

base = Path(__file__).parent
pdf_dir = base / "pdfs"
txt_dir = base / "text"
txt_dir.mkdir(exist_ok=True)

for pdf in sorted(pdf_dir.glob("*.pdf")):
    out = txt_dir / (pdf.stem + ".txt")
    try:
        reader = PdfReader(str(pdf))
        text = "\n".join((page.extract_text() or "") for page in reader.pages)
        out.write_text(text, encoding="utf-8")
        print(f"OK   {pdf.name} -> {out.name} ({len(text)} chars, {len(reader.pages)} pages)")
    except Exception as e:  # noqa: BLE001
        print(f"FAIL {pdf.name}: {e}", file=sys.stderr)
