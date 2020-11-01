#!/usr/bin/env python
  
import setuptools

setuptools.setup(
        name="runflex",
        version="0.0.1",
        author="Guillaume Monteil",
        author_email="guillaume.monteil@nateko.lu.se",
        description="runflex",
        long_description_content_type="text/markdown",
        packages=setuptools.find_packages(),
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: European Union Public Licence 1.2 (EUPL 1.2)",
            "Operating System :: OS Independent",
        ],
        python_requires='>=3.6',
        data_files=[]
)
