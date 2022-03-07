# Creating Documentation

Creating the **PPG-beats** documentation.

---

## Create and preview the documentation locally

Pre-requisites:

- **MkDocs**: _MkDocs_ is a 'static site generator that's geared towards building project documentation (from [mkdocs.org](https://www.mkdocs.org/)). Installation instructions are available [here](https://www.mkdocs.org/user-guide/installation/). _e.g._ on my Mac, I installed MkDocs using the following command: ``pip install mkdocs``

Steps:

1. **Download repository:** Download the [GitHub repository](https://github.com/peterhcharlton/ppg-beats) (which contains both the code and documentation files): Use [this link](https://github.com/peterhcharlton/ppg-beats/archive/refs/heads/main.zip) to download a ZIP file.
2. **Extract:** Extract (unzip) the ZIP file.
3. **Set current directory:** Go to Command Prompt (on Windows) or Terminal (on MacOS). Set the current directory to the extracted folder. _e.g._ on my Mac, I use: ``cd /Users/petercharlton/Documents/GitHub/ppg-beats/``
4. **Preview documentation:** Use the ``mkdocs serve`` command to view the documentation in your browser at [http://127.0.0.1:8000/](http://127.0.0.1:8000/).

## Upload the documentation to the web

Pre-requisites:

- **GitHub**: A GitHub account (and preferably some experience with GitHub).
- **Read the Docs**: A _Read the Docs_ account (although this can be easily set up, and no experience is required).

Steps:

1. **Upload to GitHub**: Commit the documentation to a GitHub repository. _e.g._ in my case, [peterhcharlton/ppg-beats](https://github.com/peterhcharlton/ppg-beats). Note that the documentation is committed to the main branch. It contains:
    - a file named _mkdocs.yml_ in the main folder
    - a subfolder named _docs_ which contains markdown documentation
    - a subfolder named _source_ which contains the source code
2. **Host on _Read the Docs_**: Host the documentation on _[Read the Docs](https://readthedocs.org/)_, allowing it to be viewed as a website. _e.g._ the _PPG-beats_ documentation is hosted by _Read the Docs_ [here](https://readthedocs.org/projects/ppg-beats/). The following are helpful for working out how to do this:
    - [Read the Docs Tutorial](https://docs.readthedocs.io/en/stable/tutorial/): Provides details of how to host content stored in a GitHub repository.
    - [MkDocs 'Deploying your docs' guide](https://mkdocs.readthedocs.io/en/stable/user-guide/deploying-your-docs/)

