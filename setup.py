from setuptools import setup, find_packages, Extension

setup(
    name='spkmeans',
    version='1.0.0',
    author= 'Talya Francus',
    description='kmeans clustering',
    install_requires=['invoke'],
    packages=find_packages(),
    ext_modules=[
        Extension(
            'spkmeans',
            ['spkmeans.c']
        )
    ]
)