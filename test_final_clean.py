#!/usr/bin/env python3
from qxub.scheduler import JobSpinner
import time

print("Testing clean spinner flow:")

print("\n1. Job submission (spinner only)")
with JobSpinner(show_message=False, quiet=False):
    time.sleep(2)

print("🚀 Job submitted successfully! Job ID: 12345")

print("\n2. Waiting with spinner on same line")
with JobSpinner("Waiting for job to start...", show_message=True, quiet=False):
    time.sleep(3)

print("Done - spinner should be cleared properly!")