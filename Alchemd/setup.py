# from setuptools import setup, find_packages

# setup(
#     name="Alchemd",
#     version="0.1.0",
#     packages=find_packages(),
#     install_requires=[
#         "openmm",  # 添加其他依赖
#     ],
#     entry_points={
#         'console_scripts': [
#             'alchemd-fep=alchemd.openmm_fep_run:main',  # 创建命令行入口点
#         ],
#     },
#     author="Your Name",
#     author_email="your.email@example.com",
#     description="A package for FEP calculations using OpenMM",
#     long_description=open("README.md").read(),
#     long_description_content_type="text/markdown",
#     url="https://github.com/yourusername/Alchemd",
#     classifiers=[
#         "Programming Language :: Python :: 3",
#         "License :: OSI Approved :: MIT License",
#         "Operating System :: OS Independent",
#     ],
#     python_requires=">=3.6",
# )
from setuptools import setup, find_packages

setup(
    name='alchemd',
    version='0.1.0',
    packages=find_packages(exclude=["example"]),
    description=(
        'An alchemical-based free energy calculation package (supporting DC-ABFE) using OpenMM'
    ),
    author="Zhe Li, Runduo Liu",
    author_email="lizhe5@mail.sysu.edu.cn",
    platforms=["linux"],
    url="https://github.com/ZheLi-Lab/Alchemd_fep_dcabfe",
    entry_points={
        'console_scripts': [
            'run_alchemd = alchemd.openmm_fep_run:main'
        ]
    },
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    package_data={
        'alchemd': ['alchemd/example']
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
