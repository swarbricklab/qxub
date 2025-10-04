#!/bin/bash
# Pylint check script for qsub_tools
cd /g/data/a56/software/qsub_tools
source venv/bin/activate
python -m pylint $(git ls-files '*.py')