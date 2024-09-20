from setuptools import setup

setup(
    name="qconda",
    version="0.1",
    py_modules=["qconda"],  # Name of the Python file without the extension
    install_requires=["click"],
    package_data={
        '': ['jobscripts/qconda.pbs'],
    },
    entry_points={
        'console_scripts': [
            'qconda = qconda:conda_submit',
        ],
    },
)

