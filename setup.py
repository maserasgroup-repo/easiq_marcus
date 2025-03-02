from setuptools import setup, find_packages

setup(
    name="pyoset",
    version="0.1",
    author="Lucía Morán-González, Albert Solé-Daura",  
    py_modules=['4points', 'get_info', 'pyosethermo', 'simpledft_marcus','utilities'],
    packages=find_packages(),
    install_requires=['pyssian @ git+https://github.com/maserasgroup-repo/pyssian.git@develop', 
                      'pyssianutils @ git+https://github.com/maserasgroup-repo/pyssian-utils.git'] # Add dependencies if needed
)
