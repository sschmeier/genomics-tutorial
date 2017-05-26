"""
Used for two things: 1. Locally for git workflow, 2. Update remote servers.
"""
from __future__ import with_statement
from fabric.api import *
from fabric.colors import *
from fabric.context_managers import cd, lcd
import os.path

HOST = '' # e.g. seb@bla.com
REMOTE_BASE_DIR = '/webapps/seb_django/www'  # absolute path, where project/repo lives
REMOTE_ERR_FILE = '/webapps/seb_django/logs/00UPDATE_203341.err'  # absolute path
REMOTE_LOG_FILE = '/webapps/seb_django/logs/00UPDATE_203341.log'  # absolute path
REPO_NAME = 'genomics-tutorial'  # basename of project
REPO_URL = 'git@github.com:sschmeier/genomics-tutorial.git'  # e.g. github url
REPO_BRANCH = 'gh-pages'  # this is the branch to clone on servers


@hosts('seb@vm010865.massey.ac.nz', 'seb@vm010944.massey.ac.nz') # only for deploy
def logs():
    """ Reading Massey server log-files and print to stdout."""
    puts(yellow("[Reading log-file]"))
    run("cat %s" % REMOTE_LOG_FILE)
    puts(yellow("[Reading err-file]"))
    run("cat %s" % REMOTE_ERR_FILE)


@hosts('seb@vm010865.massey.ac.nz', 'seb@vm010944.massey.ac.nz') # only for deploy
def deploymassey(activate_env=True, conda=None):
    """ Deploy project to Massey servers."""
    remote_dir = os.path.abspath(os.path.join(REMOTE_BASE_DIR, REPO_NAME))
    
    if activate_env:
        if conda:
            puts(yellow("[Activate conda env]"))
            run('source activate %s' %(conda))
        else:
            puts(yellow("[Activate env through ~/bin/activate]"))
            run('source ~/bin/activate')

    with settings(warn_only=True):
        if run("test -d %s" % (remote_dir)).failed:
            puts(red("[Repo %s does not exist on remote at: %s]" % (REPO_NAME, remote_dir)))
            with cd(REMOTE_BASE_DIR):
                run("git clone -b %s %s %s" % (REPO_BRANCH, REPO_URL, REPO_NAME))

    puts(yellow("[Write logs]"))
    run("echo '-----------------------------' > %s" % REMOTE_ERR_FILE)
    run("echo `date` >> %s" % REMOTE_ERR_FILE)
    run("echo '-----------------------------' >> %s" % REMOTE_ERR_FILE)
    run("echo '-----------------------------' > %s" % REMOTE_LOG_FILE)
    run("echo `date` >> %s" % REMOTE_LOG_FILE)
    run("echo '-----------------------------' >> %s" % REMOTE_LOG_FILE)

    puts(yellow("[Update repo: %s]" % REPO_NAME))
    with cd(remote_dir):
        run("git pull origin %s >> %s 2>> %s" %
            (REPO_BRANCH, REMOTE_LOG_FILE, REMOTE_ERR_FILE))

        
def git(br, to_br='master', v=None):
    """Execute local git checkout master, merge branch into master.

    Keyword arguments:
    br -- the branch that should be merged into 'to_br'
    to_br -- branch to merge to (defaults to master).
    v -- new version/tag number requested this will create a repo tag.

    Usage:
    fab git:br='new_feature',v='v1.2.5'
    """

    # co master and merge
    puts(yellow("[Checkout branch %s]"%(to_br)))
    local("git checkout %s"%(to_br))

    puts(yellow('[Merge branch "%s" into %s]'%(br,to_br)))
    local("git merge %s --no-ff" %br)

    with settings(warn_only=True):
        if v:
            puts(yellow("[Bump version: %s]"%v))
            # bump version number: project specific
            local("sed -i '' 's/v.\..\../%s/g' VERSION.txt" %v)
            # add config.json and commit
            local("git add VERSION.txt")
            local('git commit -m "Bumped to %s"' %v)

            # add tag
            puts(yellow("[Tag new version: %s]"%v))
            local('git tag -a %s'%v)


def deploy(msg, br='master'):
    """Deploy master to remotes.
    - make latexpdf
    - copy pdf to _static.
    - Safe conda packages of used env in package-list.txt
    - commit changes.
    - push master to gitlab remote.
    - push master to origin (github)
    - make clean; make html
    - commit html to gh-pages
    - push gh-pages to github/gh-pages

    Keyword arguments:
    msg -- commit message
    br -- the branch that should be pushed

    Usage:
    fab deploy:msg="This is a commit message"
    """
    # co branch
    puts(yellow("[Checkout branch %s]"%(br)))
    local("git checkout %s"%(br))

    # create new pdf
    puts(yellow('[Make latexpdf]'))
    local("make latexpdf > /tmp/latex")

    puts(yellow('[Copy pdf]'))
    local("cp ../build/latex/*.pdf _static/")

    # save package list
    puts(yellow('[Conda package list]'))
    local("conda list --export > conda-package-list.txt")

    # save conda env
    puts(yellow('[Conda env freeze]'))
    local("conda env export > conda-env-freeze.txt")

    puts(yellow('[git stage/commit changes]'))
    local("git add -u")
    local('git commit -m "%s"' %(msg))
    
    # push changes to gitlab
    puts(yellow("[Push %s to gitlab]"%(br)))
    local("git push gitlab %s"%(br))

    # push changes to github
    puts(yellow("[Push %s to github]"%(br)))
    local("git push origin %s"%(br))

    # for html to github gh-pages
    puts(yellow("[Make html for github/gh-pages]"))
    local("make clean; make html")

    puts(yellow("[Push gh-pages to github]"))
    puts(red("Will NOT add newly created content. Only already tracked content."))
    with lcd("../build/html"):
        local("git add -u")
        local('git commit -m "%s"'%(msg))
        local("git push origin gh-pages")

