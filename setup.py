from setuptools import setup, find_packages

setup(
    name='tau_project',
    version='0.1.0',
    description='A scientific simulation package for tau protein phosphorylation and aggregation',
    author='Your Name',
    author_email='your.email@example.com',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'numpy',
        'matplotlib',
        'seaborn',
        'pandas',
    ],
    tests_require=['pytest'],
    python_requires='>=3.8',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    include_package_data=True,
) 