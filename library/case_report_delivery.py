"""
What a deployment specific app did with a case report, in a shape the public Reports card can show.

A delivery is one line in the card's Delivery column - where the report was sent, how that went, and
the button to send it again. Apps answer classification's case_report_deliveries_signal with these,
so core never learns which app delivers what (the same arrangement as library/integration_status.py).
"""
from dataclasses import dataclass
from datetime import datetime
from typing import Optional


@dataclass
class CaseReportDelivery:
    """ One "Mocha: Sent, matched" row against a report on the Reports card """
    label: str  # the system the report went to
    status: str  # bootstrap contextual class - success / warning / danger / secondary
    text: str  # "Sent, matched" / "Unresolved: no Mocha extraction for 25-245-16107"
    timestamp: Optional[datetime] = None
    action_url: Optional[str] = None  # a POST that retries / resends, run by the card's action handler
    action_label: Optional[str] = None  # "Send to Mocha"
