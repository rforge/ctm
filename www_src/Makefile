## Makefile for ctm-homepage ##

## build html-files for website
build: pages-dist
	bundle exec jekyll build
	cp -r _site/* ../www/
	@echo "NOTE: To update the website, please add/commit changes in ctm/www/."

## Development only
dev: pages-dist
	bundle exec jekyll serve --config _config.yml,_config_dev.yml
	
## Rmd files in pages
pages-dist:
	make -C pages md
	
clean:
	make -C pages clean
	rm -rf _site
