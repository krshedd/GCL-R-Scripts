# GCL-R-Scripts
A repository for all GCL R Scripts. This repository is meant to store and provide version control for all shared GCL scripts. New functions should be developed and included as necessary. Older functions that previously existed in the "temp" folder will be migrated into this repository over time, so all can use them. GitHub is a collaborative tool and is meant to store not only our code, but to also allow us to keep all R script related issues and discussions in one accessible place.  
Files were initially moved off of the V: drive on 2/29/16.

## Initial cloning of repository
* Download [RStudio](https://www.rstudio.com/products/RStudio/#Desktop)
* Dowload [Git](https://git-scm.com/download/win)
* Create a [GitHub account](https://github.com/join?source=header-home)
* Enable version control in RStudio via Tools --> Global Options --> Git/SVN
* Open your git-bash.exe ("C:\\Program Files\\Git\\git.bash.exe") and follow code below, changing to your *username* and *e-mail*. Note that the "Repository URL" does not change to your username as it lives in `krshedd`'s account.
    * Configure Git with your GitHub account
        * Set your username to your GitHub username `git config --global user.name "krshedd"`
        * Set your e-mail address to what you use for GitHub `git config --global user.email "kyle.shedd@alaska.gov"`
        * Verify changes `git config --global --list`
    * Set up SSH key for RStudio to communicate with GitHub
        * Change directory to where you want to save key `cd /c/Users/krshedd`
        * Make directory to save key in `mkdir ".ssh"`
        * Create key `ssh-keygen -t rsa -b 4096 -C "kyle.shedd@alaska.gov"`
        * Enter file in which to save the key `/c/Users/krshedd/.ssh/id_rsa`
        * Create passphrase for key
        * Re-enter passphrase for key
        * Open ssh-agent `eval "$(ssh-agent -s)"`
        * Add key to ssh-agent `ssh-add /c/Users/krshedd/.ssh/id_rsa`
        * Enter passphrase for key
        * Copy SSH key `clip < /c/Users/krshedd/.ssh/id_rsa.pub`
        * Follow [instructions](https://help.github.com/articles/adding-a-new-ssh-key-to-your-github-account/) to add SSH key to your GitHub account online
        * Test your SSH connection `ssh -T git@github.com`
        * Are you sure you want to continue connecting `yes`
        * After confirmation, close git.bash.exe
* Initial clone of `GCL-R-Scripts` repository on GitHub to your local machine
    * Open a new project in RStudio via File --> New Project... --> Version Control --> Git
        * Repository URL: <https://github.com/krshedd/GCL-R-Scripts.git>
        * Project directory name autofills
        * Create project as subdirectory of: wherever you want on your C: drive (i.e. "C:\\Users\\krshedd\\Documents\\R")
        * Close project; you now have a local copy of all files from the `GCL-R-Scripts` repository.
* As a final step, you need a `Functions.GCL.R` R-script that will source all of the R-scripts from your local directory. Note that you will need to change `directory` to the path to your local "cloned" repository.

<pre><code>#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# This function sources GCL functions from my C drive (GitHub clone)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
directory <- "C:/Users/krshedd/Documents/R/GCL-R-Scripts"
functions <- list.files(path = directory, pattern = "\\.R$", full.names = TRUE, recursive = FALSE, ignore.case = TRUE)  # read in all .R scripts
EMPTY <- sapply(functions, source)
rm(directory, functions, EMPTY)
print(as.vector(lsf.str()))
</code></pre>
* Once you have created the above `Functions.GCL.R`, you will source to it directly to access all `GCL-R-Scripts`  
`source("C:/Users/krshedd/Documents/R/Functions.GCL.R")`

## Normal usage
Best practice is to **pull** frequently (i.e. beginning of each day).

* Open your `GCL-R-Scripts.Rproj` in RStudio
* In the upper right panel select **Git**
* Confirm that you are on the `master` branch (upper right)
* Select the blue down arrow for **Pull**. A window will indicate whether you are `already-up-to-date` or inform you of changes since your last pull
* Close `GCL-R-Scripts.Rproj`
* You now have the most up to date version of all R scripts on your local C: drive

## Releases
Periodically, after large changes to the `GCL-R-Scripts` repository, a specific commit will be tagged as a *release*. Each release will have a version number and a thorough set of release notes indicating:  
1. New Functions  
2. Function Updates  
3. Bug fixes  
Release notes serve to update you periodically of major changes so you don't have to go through each commit specifically to note changes.

## Branches
The `GCL-R-Scripts` repository has two main branches: `master` and `develop`. The thought is to loosely follow this [successful Git branching model](http://nvie.com/posts/a-successful-git-branching-model/). The `master` branch is what we should all be pulling from and sourcing to. Ideally, the `develop` branch is where new functions are added and changes are made to existing scripts between *releases*. When it is time for a new *release* the `develop` branch is merged back to the `master` branch once code has been vetted and tested. If there are bugs in a *release* that need immediate attention, a temporary `hotfix` branch should be created from the `master` to address the bug.

## Issues
One of the most powerful features of GitHub is tracking *Issues*. When someone has a problem with one of our functions, finds a bug, or wants a new feature, it is best to create a *Issue*. This allows us to:  
1. Document the problem a user had with a time-stamp,  
2. Discuss potential solutions,  
3. Come to a consensus solution,  
4. Assign a person to resolve the issue, and  
5. Write code that solves the problem and closes the *Issue*  
This is a much better solution than e-mailing Jim or individual users. This allows everyone to know what the problem is, who had it, and if it has been fixed. *Issues* can be [closed by commits](https://help.github.com/articles/closing-issues-via-commit-messages/), and should not be closed until code is written to solve the *Issue*.

## Sourcing scripts from Remote Desktop
If you need to source R scripts, but do not have access to your local C: drive (i.e. when you are on the remote desktop servers), you can source scripts individually by their unique https url. All R scripts can be accessed by their raw https url which follows the formula "https://raw.githubusercontent.com/*username*/*repository*/*branch*/*filename*", for example <https://raw.githubusercontent.com/krshedd/GCL-R-Scripts/master/AttributesToIDs.GCL.r>.  
Another way to find the raw https url is to navigate within the repository on GitHub to the function, select *Raw*, and copy the url.  
Once you have the url for the function you want to source there are two options to source it. The first one is easiest (less code), but the second one will always work (requires the package "RCurl").  

<pre><code>source(https://raw.githubusercontent.com/krshedd/GCL-R-Scripts/master/AttributesToIDs.GCL.r)`</code></pre>

<pre><code>while(!require(RCurl)){install.packages("RCurl")}
eval(expr = parse(text = getURL("https://raw.githubusercontent.com/krshedd/GCL-R-Scripts/master/AttributesToIDs.GCL.r", ssl.verifypeer=FALSE) ))</code></pre>

## QC
As of release [QC For Everyone v1.4.0](https://github.com/krshedd/GCL-R-Scripts/releases/tag/v1.4.0) the standard QC R script run by project leaders will work for both SNPs and uSATs and is available on the `master` branch. This is a separate script all of its own, not a full-blown function. Weird things always seem to happen during QC, so rather than have a black box function, it was decided to create a dedicated script that project leaders would then run line by line so that they can pinpoint problems if/when they arise. The idea is that when you need to perform a QC, you will do a fresh `pull` from the `master` branch in order to have the latest version of `QC.R` in your remote repository (i.e. `C:\Users\krshedd\Documents\R\GCL-R-Scripts`). You will then **copy/paste** this file to the QC project folder (i.e. `V:\Lab\Genotyping\SNP Projects\Sockeye\Project S161 Chignik Inseason 2016\QC`), modify the **Arguments**, and run the script from there. You should not have to modify anything past **GO!**, but if the function bombs on you, you will at least know where you are. This arrangement will also allow you to save your project specific QC script on the `V drive` when you are done.