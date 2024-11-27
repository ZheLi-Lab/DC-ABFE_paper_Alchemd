from setuptools import setup, find_packages

setup(
    name='alchemd_core',
    version='0.1.0',
    packages=find_packages(exclude=["example"]),
    description=(
        'An alchemical-based free energy calculation package (supporting DC-ABFE) using OpenMM'
    ),
    author="Zhe Li, Runduo Liu",
    author_email="lizhe5@mail.sysu.edu.cn",
    platforms=["linux"],
    url="https://github.com/ZheLi-Lab/DC-ABFE_paper_Alchemd",
    entry_points={
        'console_scripts': [
            'run_alchemd = alchemd_core.openmm_fep_run:main'
        ]
    },
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    package_data={
        'alchemd_core': ['alchemd/example']
    },
    python_requires = '>=3.9,<3.11',
    install_requires=[
        'openmm==7.7',
        'openmmtools==0.21.5',
        'tqdm',
        'alive-progress',
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: Other/Proprietary License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
        "Intended Audience :: Science/Research",
    ]
)
