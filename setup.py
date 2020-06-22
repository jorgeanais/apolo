from setuptools import setup, find_packages

setup(
    name='apolo',
    version='0.2',
    description='Clustering utilities for astronomical datasets',
    author='Jorge Anais',
    url='https://github.com/jorgeanais/apolo',
    author_email='jrganais@gmail.com',
    license='MIT',
    packages=find_packages(),
    install_requires=['astroquery', 'numpy', 'matplotlib', 'astropy', 'scipy', 'seaborn', 'hdbscan', 'scikit-learn'],
    package_data={
        '': ['*.md'],
    },
    zip_safe=False,
)