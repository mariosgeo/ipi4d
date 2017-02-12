========================================================================

This readme file provides guidelines and good practices for the manage-
ment of the IP4DI repository using Git.

Francois Lavoue, 23 November 2016
Last updated 12 Feb. 2017

========================================================================


=======================
Outline of this readme:

- Organisation and management of the IP4DI repository
  - Branches of the IP4DI repository
  - How to work on the repo

- Generalities about Git     (Git beginners may start here)
  - Useful links about Git
  - Git in a nutshell
  - Important terms when speaking about Git
  - Architecture of a Git repository
  - Some things that are good to know
  - Good practices
  - Bad practices


===================================================
Organisation and management of the IP4DI repository

-----
Branches of the IP4DI repository

 - master       # common branch, can be updated by 'master' members only
 - stable       # should be identical to master, can be updated only by
                  the owner for public distribution on GitHub
 - dev_IGI_FL   # development branch associated to Image-Guided Inver-
                  sion (IGI), and to the "no-GUI" version of the code,
                  can be updated by FL only.
 - appli_Solfatara_FL   # application branch dedicated to the applica-
                          tion of the code to the Solfatara dataset.

Development branches may -or may not- be intended to be merged to the
master branch, or to any other dev. branch in the long term. Within a
group, each member should have its own dev. branch, or several ones
associated to different developments. A typical name for a development
branch should be something like:
      dev_$short-descriptive-name_$DEVELOPER-INITIALS
ex:   dev_IGI_FL

Merges to the master branch should only concern debugging or developments
that are thought to be useful for all users. They are only performed by
'master' members of the project (people who have enough experience with
this particular code to not mess things up).

Application branches can be useful to keep track of fine tuning of the
code for specific applications. They can also be a good tool in the view
of publishing reproducible research.

NB: while the stable branch is made public on GitHub, the development
    branches are hosted on GitLab.com which, unlike GitHub, enables the
    management of private repositories for free.


-----
How to work on the repo

In practice, people should the following steps:

1) Create an account on www.gitlab.com. You may also register your SSH
   keys (see section "good to know" below).

2) Ask the owner of the repository (currently Francois Lavoue) to give
   you access to the repo.

3) Start by cloning the repository on your own machine:
     cd YOUR_WORKING_DIR
     git clone git@gitlab.com:lavouef/IP4DI.git .

4) Then THE FIRST THING TO DO BEFORE MAKING ANY CHANGE is to create your
   own branch:
     git checkout -b dev_your-dev-topic_your-initials
   and push it to the remote repository on gitlab.com (or not, if you
   prefer keeping it hidden from the other group members):
     push -u origin dev_your-dev-topic_your-initials

5) You can also create new projects (see https://gitlab.com).



=========================
Generalities about Git

-----
Users new to Git may refer to these useful links:

http://git-scm.com/docs/gittutorial      (to begin with Git)
http://git-scm.com/docs/gitcvs-migration (the centralized way of using Git, as we do here)
http://git-scm.com/docs/gitworkflows     (good practices)
http://git-scm.com/docs/giteveryday      (useful commands according to user level)
http://git-scm.com/docs/gitglossary      (definition of Git terms)
http://betterexplained.com/articles/aha-moments-when-learning-git/  (Git branches are VIRTUAL directories)


-----
Put in a nutshell, Git is a version control system (VCS) that enables

- to work in a collaborative and organized way on source codes (or any
  type of document) by keeping track of the changes (versioning).
- to backup all previous versions of the code, such that we can go back
  to an earlier version in case of bug.
- to manage projects and project members online (www.gitlab.com).


-----
Important terms when speaking about Git are

- 'repositories' which contain the source code files and a '.git' folder
  that enables the management of these files (in particular, it contains
  the whole history of the source files).
- within these repositories, 'branches' are variants of the code develo-
  ped in parallel by different people, or by the same developer for dif-
  ferent aims.
- a 'version' is the state of a given branch at a given time.


-----
Architecture of a Git repository

  The default branch in a Git tree is generally the 'master' branch.
A 'stable' branch may also exist and contain the version of the code that
should be distributed to non-Git users and used for computations.
  From the master branch, we can derive development branches where modi-
fications are implemented and tested. In principles, these development
branches are thought to be merged back to the master branch once develop-
ments are finished.


-----
Some things that are good to know:

- Git branches do NOT coincide with physical folders! This is a bit
  weird at first sight but you'll get used to it...

- fetch vs. pull:
  - git fetch is the command that says
      "bring my local copy of the remote repository up to date"
    (it does not changes the code you are working on, but just your
    local copy of the remote repository, with which you compare your
    code when doing "git status")
  - git pull says
      "bring changes from the remote repository into my own code"
    (effectively changes the code you are working on)
  (for details, see http://stackoverflow.com/questions/292357/what-are-the-differences-between-git-pull-and-git-fetch)

- If you are new on git.dias, you may register the SSH key(s) of your
  machine(s) on git.dias (under Profile Settings > SSH Keys), such that
  git.dias does not ask you your login and password at each pull or push
  push command. Then make sure to use the SSH protocol to clone, pull or
  push. "git remote -v" should display a remote url of the form
    git@git.dias.ie:$GIT_GROUP/name-of-repo
  If you have instead
    https://git.dias.ie/$GIT_GROUP/name-of-repo
  then change the remote url using
    git remote set-url origin git@git.dias.ie:$GIT_GROUP/name-of-repo


-----
Good practices:

- BE CAREFUL to which branch you are working on (ie work on YOUR branch).
  As Git branches are VIRTUAL and can be switched from one to the other
  in the same physical directory, it is quite easy to mess up branches...
  "git branch" or "git status" with tell you on which branch you are.

- Even as a 'master' member, ALWAYS CONSULT the other $GIT_GROUP members
  BEFORE MERGING anything in the MASTER branch.

- Please provide explicit commit messages (not just "commit" or "add good
  stuff"). Commit messages should have the following structure:

    === Example of a commit to the IP4DI package ===

    Major update of the no-GUI master branch.                           # short descriptive title)

    Main modifications since last master update include:                # longer description and details if needed...

    - no-GUI time-lapse inversion
    - consideration of topography in the no-GUI version
    - further developments of image-guided inversion (still in progress...)
    - development of plotting tools

    Francois Lavoue, 22 Nov. 2016                                       # Author Name, Date

    -----------------------------                                       # separate description from list of modified files

      On branch master                                                  # branch name and list of modified files
      Changes to be committed:                                          # (you just have to uncomment some lines in the
        modified:   00README.txt                                        #  default commit message to make them appear)
        modified:   ip4di_2d_no-GUI_IGI/00README.txt                    # (do not include untracked files, or files that
        deleted:    ip4di_2d_no-GUI_IGI/0readme_git.txt                 #  are not part of this particular commit)
        modified:   ip4di_2d_no-GUI_IGI/forward_parameters.m
        ...

    === End example ===

- Commit your changes regularly (and frequently), in order to (1) backup
  your work, (2) be able to come back to a recent bug-free version in
  case of bugs, and (3) avoid conflicts during future merges of a parent
  branch (eg master into your dev. branch).

- In the same order of idea, one commit should be associated with only
  one "change". This does not mean that you have to commit each time you
  change a comma: this "change" may concern several source files but all
  of these changes should target the same aim.
  Examples: (1) change a scalar variable from single to double precision
  in all routines where it appears, or (2) make an array allocatable in
  the routine where it is defined first, but not (1+2) at the same time:
  these two changes should be the subject of two different commits.
    Do this by choosing the files to be committed with "git add file1
  file2" rather than "git add .".
    This makes the cherry-picking of the different changes much easier.

- By default, Git gives weird SHA-1 names to commits. It is good practice
  to give more understandable names (e.g. vX.X) to important commits in
  order to be able to retrieve them easily in the logs later on.
  Do this through the "git tag" command.

- Important: Try to stick to the original code as much as possible in
  your developments, in terms of inputs/outputs, variable names, and
  "spirit", such as to ensure back-compatibility and avoid conflicts
  during merges and cherry-picking.


-----
Bad practices:

- As said earlier, DO NOT MERGE ANYTHING INTO MASTER WITHOUT CONSULTING
  the other $GIT_GROUP members...

- AVOID going back to a previous version using "git reset --hard": do
  NOT do this if you are sharing your branch with other people who have
  copies of the old commits, because using a hard reset will force them
  to resynchronize their work with the newly reset branch. For a solu-
  tion that explains in detail how to safely revert commits without
  losing work with a hard reset, see
http://stackoverflow.com/questions/4114095/revert-git-repo-to-a-previous-commit/4114122#4114122.

    Ultimately, if you messed up things and cannot go back properly to
  erase the mess from Git history, simply correct your mistake and com-
  mit a debugged version with an explicit message specifying the issue
  (and the fact that the previous commit is bugged). We better have an
  extra commit and a dirty but consistent history than a faked, broken,
  history that does not enable to restore earlier versions.

- New branches should not be created from stable, but from master, or
  from other dev. branches, such that stable only interacts with master,
  and only to merge changes from master to stable. This is both to pro-
  tect the stable branch from possible mistakes and to ensure you get
  the latest master version in case it would not have been merged with
  stable yet.

