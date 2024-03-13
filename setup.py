from setuptools import setup, find_packages

from pkg_resources import parse_requirements

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

with open("requirements.txt", encoding="utf-8") as fp:
    install_requires = [str(requirement) for requirement in parse_requirements(fp)]

setup(
    name="PyArc",
    version="1.0.0",
    author="Siyuan Xu",
    author_email="2016302540149@whu.edu.cn",
    description="A python Package for computing absorption"
                " coefficients and radiative recombination rates in semiconductors",
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT License",
    url="https://github.com/xsy1999/PyArc/blob/main/README.md",

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],

    packages=find_packages(),
    install_requires=install_requires,
    python_requires='>=3.6',
)