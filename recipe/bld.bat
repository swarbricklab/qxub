@echo on

:: Install the package
%PYTHON% -m pip install . -vv
if errorlevel 1 exit 1
