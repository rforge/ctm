################################################################################

# WEBSITE: Transformation Models

################################################################################

Source-folder for building the html-files with Jekyll 
Jekyll: https://jekyllrb.com
Theme: https://phlow.github.io/feeling-responsive/

Development only: jekyll serve --config _config.yml,_config_dev.yml
for more information see: _config_dev.yml

# WORK FLOW
_config.yml
	adapt configuration data
	(_config_dev.yml for development only)

pages/*
	change Rmd/md-files
	make md

pages/pages-root-folder/*
	index.md: Provided that the file has a front matter section


_layouts/: templates for frontpage/page/post
	# {{content}}-float used to inject content from md-files

_includes/: partials that can be mixed/matched '_layouts' and 'posts' to facilitate reuse

_data/: contains yml-files for variables, such as authors, etc ... also navigation



_drafts/: unpublished posts (format: title.MARKUP)


_sass/
	partials imported into "main.scss" => processed into stylesheet "main.css" (defines style of site)

_posts/
	dynamic content (format: YEAR-MONTH-DAY-title.MARKUP)

_site/	# generated sites from jekyll
	cp _site/* ../www/

