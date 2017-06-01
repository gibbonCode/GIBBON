# Project Pages

Project Pages is Jekyll Template specifically geared towards collaborative science. For more information, click [here](https://github.com/projectpages/project-pages/wiki/).

# Nav Bar Jumbles

If you have seemingly random pages popping up on your Nav Bar recently, this is due to the fact that GitHub/Jekyll changed a fundemental rule they used to render pages. 

## CAUSE:
It used to be that if a markdown file didn't have `---` frontmatter at the beginning, it wasn't rendered as a page. This was changed very recently (like in the last 2 days) so that every markdown file anywhere no matter what gets rendered as a page.  

## FIX:

1) Go to:

`project-pages/plugin/projector/` or `yourreponame/plugin/projector/` and delete the `README.md` file. This can be done graphically for the non-Git-savvy by simply going to your:

GitHub account -> Your Profile -> Repositories -> Project-Pages/Your Repo -> Plugin -> projector 

and clicking on the files, then clicking on the "thrash can / delete this file" icon on the top right corner of the file.

2) Go to:

`project-pages/css/theme/` or `yourreponame/css/theme/` and delete the `README.md` file. This can be done graphically for the non-Git-savvy by simply going to your:

GitHub account -> Your Profile -> Repositories -> Project-Pages/Your Repo -> Plugin -> projector 

and clicking on the files, then clicking on the "thrash can / delete this file" icon on the top right corner of the file.
