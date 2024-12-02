from setuptools import setup, find_packages
version = "0.0.1"
setup(
    name="moldscript",
    packages=find_packages(exclude=["tests"]),
    package_data={"molDscript": ["templates/*"]},
    version=version,
    license="MIT",
    description="An Automated Workflow for Quantum Mechanical Derived Descriptors: MOLDSCRIPT",
    long_description="Documentation in Read The Docs: XXX",
    long_description_content_type="text/markdown",
    author="Shree Sowndarya S. V., Jake King",
    author_email="svss@colostate.edu, jake.king@colostate.edu",
    keywords=[
        "molecular descriptors",
        "electronic parameters",
        "steric parameters",
        "workflow",
        "computational chemistry",
        "cheminformatics",
        "quantum mechanics",
        "DFT",
        "automated",
    ],
    url="https://github.com/patonlab/molDscript",
    download_url=f"https://github.com/patonlab/molDscript/archive/refs/tags/{version}.tar.gz",
    classifiers=[
        "Development Status :: 3 - Alpha",  # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        "Intended Audience :: Developers",  # Define that your audience are developers
        "Topic :: Software Development :: Build Tools",
        "License :: OSI Approved :: MIT License",

    ],
    install_requires=[
        "pandas>=2.0.2",
        "cclib @ git+https://github.com/cclib/cclib.git",
        "dbstep",
        "openbabel",
        'rdkit',

    ],
    python_requires=">=3.0",
    include_package_data=True,
)
