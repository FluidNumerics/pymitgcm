import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="pymitgcm-fluidnumerics_joe", # Replace with your own username
    version="0.0.1",
    author="Joe Schoonover",
    author_email="joe@fluidnumerics.com",
    description="A package for post-processing MITgcm data in python",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/FluidNumerics/pymitgcm",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache-2.0 License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
