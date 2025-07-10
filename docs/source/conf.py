# -*- coding: utf-8 -*-
#
# This file is execfile()d with the current directory set to its containing dir.
#
# The contents of this file are pickled, so don't put values in the namespace
# that aren't pickleable (module imports are okay, they're removed automatically).
#
# All configuration values have a default value; values that are commented out
# serve to show the default value.

#import sys, os

# If your extensions are in another directory, add it here. If the directory
# is relative to the documentation root, use os.path.abspath to make it
# absolute, like shown here.
#sys.path.append(os.path.abspath("some/directory"))

# General configuration
# ---------------------

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named "sphinx.ext.*") or your custom ones.
extensions = [
  "sphinx.ext.autodoc",
  "sphinx.ext.viewcode",
  "sphinx.ext.intersphinx",
  "sphinx_copybutton",
  "sphinx_design"
]

# Add any paths that contain templates here, relative to this directory.
templates_path = ["_templates"]

# The suffix(es) of source filenames.
# You can specify multiple suffix as a list of string:
# source_suffix = ['.rst', '.md']
source_suffix = ".rst"

# Be strict about any broken references:
nitpicky = True

'''
see this: https://techwriter.documatt.com/2021/sphinx-conf-py-tips.html
ignore allows us to avoid things like:
apipy.rst:20: WARNING: py:class reference target not found: numpy array
'''
nitpick_ignore = [('py:class', 'numpy array')]

# The encoding of source files.
#source_encoding = "utf-8-sig"

# The master toctree document.
master_doc = "index"

# General substitutions.
project = "pressio-demoapps"
copyright = u"2021, National Technology & Engineering Solutions of Sandia, LLC (NTESS)"


# The default replacements for |version| and |release|, also used in various
# other places throughout the built documents.
def get_version():
  with open("../../version.txt") as version_file:
    return version_file.read().strip()

# The full version, including alpha/beta/rc tags.
release = get_version()

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
today_fmt = "%B %d, %Y"

# List of documents that shouldn't be included in the build.
#unused_docs = []

# List of directories, relative to source directories, that shouldn't be searched
# for source files.
#exclude_dirs = []

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = "tango"
pygments_dark_style = "monokai"



# Options for HTML output
# -----------------------

html_theme = "furo"

html_theme_options = {
  "sidebar_hide_name": False,
  "light_css_variables": {
    "color-brand-primary": "#336790",  # blue
    "color-brand-content": "#336790"
  },
  "dark_css_variables": {
    #"color-brand-primary": "#E5B62F"  # yellow
    #"color-brand-content": "#E5B62F",
    "color-brand-primary": "#F39C12",  # orange
    "color-brand-content": "#F39C12",
  },
}


html_sidebars = {}

# The style sheet to use for HTML and HTML Help pages. A file of that name
# must exist either in Sphinx' static/ path, or in one of the custom paths
# given in html_static_path.
#html_style = "default.css"

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = project + "\n" + "v"+release

# The name of an image file (within the static path) to place at the top of
# the sidebar.
#html_logo = None

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
#html_static_path = ['_static']

# If not '', a "Last updated on:" timestamp is inserted at every page bottom,
# using the given strftime format.
html_last_updated_fmt = "%b %d, %Y"

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
#html_use_smartypants = True

# Custom sidebar templates, maps document names to template names.
#html_sidebars = {}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_use_modindex = True

# If true, the reST sources are included in the HTML build as _sources/<name>.
#html_copy_source = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# If nonempty, this is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = ''

# Output file base name for HTML help builder.
htmlhelp_basename = "pressiodemoappsdoc"


# Options for LaTeX output
# ------------------------

# The paper size ("letter" or "a4").
#latex_paper_size = "letter"

# The font size ("10pt", "11pt" or "12pt").
#latex_font_size = "10pt"

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, document class [howto/manual]).
latex_documents = [
        ("index", "pressio-demoapps.tex", "pressio-demoapps Doc", "Francesco Rizzi", "manual"),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
#latex_logo = None

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# Additional stuff for the LaTeX preamble.
#latex_preamble = ''

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_use_modindex = True
