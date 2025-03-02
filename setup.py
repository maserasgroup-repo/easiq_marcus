from setuptools import setup, find_packages

setup(
    name="easiq_marcus",
    version="0.1",
    author="Lucía Morán-González, Albert Solé-Daura, Feliu Maseras",  
    py_modules=['4points', 'get_info', 'get_results','utilities'],
    packages=find_packages(),
    install_requires=['pyssian @ git+https://github.com/maserasgroup-repo/pyssian.git@develop', 
                      'pyssianutils @ git+https://github.com/maserasgroup-repo/pyssian-utils.git'] # Add dependencies if needed
)
