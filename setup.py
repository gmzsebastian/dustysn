from setuptools import setup

setup(
    name='dustysn',
    version='0.0.2',
    author='Sebastian Gomez',
    author_email='sebastian.gomez@austin.utexas.edu',
    description='Package to fit the SEDs of dusty supernovae.',
    url='https://github.com/gmzsebastian/dustysn',
    license='MIT',
    python_requires='>=3.6',
    packages=['dustysn'],
    license_files=["LICENSE"],
    include_package_data=True,
    package_data={'dustysn': ['ref_data/**/*', 'ref_data/*']},
    install_requires=[
        'numpy',
        'matplotlib',
        'astropy',
        'scipy',
        'emcee'
    ]
)
