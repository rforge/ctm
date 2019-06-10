################################################################################

# WEBSITE: Transformation Models
Source-folder for building the html-files with Jekyll 

################################################################################

# TODOs: 
- News: move widgets
- People: Gallery, more dynamic

# --------------- WORK FLOW  ---------------
 change .md / .Rmd files in pages, /pages-root-folder to change content

 include images in images/, either specify in .md as yml-variable (only title needed)
 or include manually by: {{ site.urlimg }}imagetitle

 pages/
	make md

 jekyll serve

 cp -r _site/* ../www/

# ------------------------------------------ 

# DEVELOPMENT ONLY (more infos: _config_dev.yml)

$ jekyll serve --config _config.yml,_config_dev.yml

# FILES
_config.yml
	configuration data
	(_config_dev.yml for development only)

_data/
	navigation.yml: website outline

pages/*
	change Rmd/md-files
	make md

pages/pages-root-folder/*
	index.md: content and figure-setting of the front-page

_layouts/: templates for frontpage/page/post
	# {{content}}-float used to inject content from md-files

_includes/: partials that can be mixed/matched '_layouts' and 'posts' to facilitate reuse

_drafts/: unpublished posts (format: title.MARKUP)

_sass/
	partials imported into "main.scss" => processed into stylesheet "main.css" (defines style of site)

_site/	# generated sites from jekyll
	cp _site/* ../www/

# BLOG
_config.yml:	defines blogurl, pagination, posturls
_posts/
	dynamic content of (blogurl): format: YEAR-MONTH-DAY-title.MARKUP
_includes/
	listing of blog-content

# FURTHER INFO
Jekyll: https://jekyllrb.com
Theme: https://phlow.github.io/feeling-responsive/
