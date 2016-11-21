Sphinx course doc
===========

## GITHUB

### Setup directory structure

```bash
$ mkdir myproject
$ cd myproject

# directory for source
$ mkdir docs
# directory for build
$ mkdir html

$ tree -L1
.
├── docs
└── html
```

Set up git in source dir:

```bash
$ cd docs
$ git init

# Add git remote
$ git remote add git@github.com:username/myproject.git
$ touch README.md
$ git add README.md
$ git commit -m "Init."
$ git push origin master

# Init Sphinx project
$ sphinx-quickstart
```

Set up the build directory:

```bash
$ cd ..
$ git clone git@github.com:username/myproject.git html
$ cd html
$ git symbolic-ref HEAD refs/heads/gh-pages
$ rm .git/index
$ git clean -fdx
```

### Makefile changes

We need to make some changes in the `Makefile` produced by `sphinx-quickstart`.
Change the build directory to one level up, which will create a html directory
upon a `make html` there automatically.

```bash
BUILDDIR      = ..
```

Change the `Makefile` clean directive to the following, which will just delete
the content in the build directory.

```bash
.PHONY: clean
clean:
	rm -rf $(BUILDDIR)/html/*
    ```

### The workflow

Make changes in the source directory `docs`. Commit changes made here to the
`master`-branch and commit freely. Push to remote `origin/master` if required.

```bash
$ git add index.rst
$ git commit -m "An update"
$ git push master origin
```

Once ready to rebuild the html, run `make html`. Change into the build directory
`../html` and commit changes to `gh-pages`-branch. Push changes to remtoe `origin/gh-pages`-branch.

```bash
$ make html
$ cd ../html
$ git add -A
$ git commit -m "New html-build"
$ git push origin gh-pages
```

Done!

## GITLAB

We also set up this repo in gitlab to make use of their pages as well. Create file `.gitlab-ci.yml` with the following content in the source directory `docs`:

```bash
image: alpine

pages:
  script:
  - apk --no-cache add py-pip python-dev
  - pip install sphinx
  - pip install sphinx_rtd_theme
  - apk --no-cache add make
  - make html
  - mv ../html/ public/
  artifacts:
    paths:
    - public
  only:
  - master
```

Then we add the gitlab remote:

```bash
git remote add gitlab git@gitlab.com:username/myproject.git
git push gitlab master
```

Done!


