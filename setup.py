#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages


requirements = ['snakehelp',
                'numpy',
                'pandas',
                'kaleido',
                'numpy_alignments',
                'snakehelp',
                'graph_read_simulator>=0.0.19',
                'tabulate'
                ]

test_requirements = ['pytest>=3', "hypothesis"]

setup(
    author="Ivar Grytten",
    author_email='ivargry@ifi.uio.no',
    python_requires='>=3.10',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.8',
    ],
    description="..",
    install_requires=requirements,
    license="MIT license",
    include_package_data=True,
    keywords='mapping_benchmarking',
    name='mapping_benchmarking',
    packages=find_packages(include=['mapping_benchmarking', 'mapping_benchmarking.*']),
    test_suite='tests',
    tests_require=test_requirements,
    version='0.0.1',
    zip_safe=False,
)
