#!/bin/bash
cd /g/data/a56/software/qsub_tools
source venv/bin/activate
echo "Line 51:"
sed -n '51p' qxub/config_manager.py
echo "Line 119:"
sed -n '119p' qxub/config_manager.py
echo "Line 121:"
sed -n '121p' qxub/config_manager.py
echo "Line 272:"
sed -n '272p' qxub/config_manager.py
