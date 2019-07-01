################################################################################

# WEBSITE: Transformation Models
Source-folder for building the html-files with Jekyll 

################################################################################

# TODOs: 
- People: Gallery, more dynamic

# --------------- WORK FLOW  ---------------

 change content of page in pages/*
	in md / Rmd files in pages/* OR pages/pages-root-folder/* for front-page
 
 for changes in pages/Rmd files:
	$ cd pages
	$ make md 

 include images in images/* 
	either specify image-title as yml-variable (in html/md-file) OR include in html/md-file as: {{ site.urlimg }}imagetitle
	# Example: <img src="{{ site.urlimg }}img.title"/>
	# Markdown-Example for posts ![Image Text]({{ site.urlimg }}image.jpg)

 include documents in docs/* 
	either specify image-title as yml-variable (in html/md-file) OR include in html/md-file as: {{ site.urlimg }}imagetitle
	# Example: <img src="{{ site.urlimg }}img.title"/>
	# Markdown-Example for posts ![Image Text]({{ site.urlimg }}image.jpg)

 new entry in "News": 
	add new file _posts/YEAR-MONTH-DAY-title.md

 build html-files:
 	$ bundle exec jekyll build
	$ cp -r _site/* ../www/


# FOR DEVELOPMENT ONTLY: $ bundle exec jekyll serve --config _config.yml,_config_dev.yml (more infos in _config_dev.yml)


# --------------- FILE-STRUCTURE ---------------

_config.yml
	configuration data: 
	(_config_dev.yml for development only)

_data/
	navigation.yml: website outline
	authors.yml: Feed for people.md
	
pages/*
	Rmd/md-files

pages/pages-root-folder/*
	index.md: content and figure-setting of the front-page

images/
	image-folder

docs/
	document-folder

_layouts/
	templates for frontpage/page/post
	COMMENT: {{content}}-float used to inject content from md-files

_includes/
	two types of partials: "_includes" are used for templates, "includes" can be used in pages or posts

_drafts/:
	unpublished posts (format: title.MARKUP)

_sass/
	partials imported into "main.scss" => processed into stylesheet "main.css" (defines style of site)

aseets/
	css/ used to format content of webpage
	img/ not changing images of webpage e.g logo

_site/	
	generated sites from jekyll


# --------------- SETUP: NEWS-PAGE (Blog-style) ---------------

_config.yml
	defines blogurl, pagination, posturls

_posts/
	dynamic content of (blogurl): format: YEAR-MONTH-DAY-title.MARKUP

_posts/_drafts
	unpublished content

_includes/
	listing of blog-content
	_meta_information: link next/previous post

change of naming:
	change _config.yml: blogurl, _data/navigation: link, 
	_data/language (for post-buttons)


# --------------- FURTHER INFO ---------------
Jekyll: https://jekyllrb.com
Theme: https://phlow.github.io/feeling-responsive/
