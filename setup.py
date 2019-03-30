# Setup.py python distribution script to install Macrocomplex-Builder in the computer

from distutils.core import setup

setup(
   name='4SMacroBuilder',
   version='0.0.1',
   description='MacrocomplexBuilder is a python program designed to 	generate macrocomplex structures from simple pair interactions',
   long_description=open('README.md').read(),
   author='B. Pau, C. Altair, S. Natàlia',
   url='https://github.com/Altairch95/4SMacroBuilder',
   packages=['4SMacroBuilder'],  
   install_requires=['biopython >= 1.73.0',
		                 'numpy >= 1.16.1'], 
   license='LICENSE.txt',
   classifiers=[
                "Programming Language :: Python :: 3",
	              "License :: OSI Approved :: MIT License",
                "Operating System :: OS Independent"],
    scripts=[
              4SMacroBuilder/MBlauncher.py ],
)
