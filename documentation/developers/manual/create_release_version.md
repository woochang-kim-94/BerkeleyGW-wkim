# BerkeleyGW: creating a release version

To create a release version of BerkeleyGW, follow the steps below:


1. Make sure you have a clean working directory of BerkeleyGW that you would
   like to make a release version from. Make sure you are on the `master` branch, and
   that you are up-do-date with your `origin` and `upstream` repos (i.e., do a `git sync`)

2. Go to the root folder of BerkeleyGW. Make sure that `Common/version_base.h`
   has the correct version number that you would like to release. You may
   likely need to increase the version number. You don't need to commit any
   change at this point.

3. Go to `documentation/user_manual/website/mkdocs.yml` and make sure that the version
   string for the manual (e.g., "BerkeleyGW 3.0") is correct.

4. Take a look at the script `dev-scripts/create_release_version.sh` to see
   what may need to be customized. In particular, make sure that the `VERSION`
   variable in the script is correct and matching that in
   `Common/version_base.h`. You may also need to change the hard-coded
   commands, such as `SEDEXE` and `GREPEXE`. Importantly, also check that all
   the folders that need to be removed from the release version of BerkeleyGW
   are properly included in the script. Again, you don't need to commit any
   change at this point.

5. After the version file and the script was corrected, run
   `./dev-scripts/create_release_version.sh`. In case of success, a file
   `BerkeleyGW-${VERSION}.tar.gz` should have been created. If not, there
   are likely missing `*_INTERNAL_ONLY` tag.

6. Extract the tarball into a clean directory. Make sure you can make the code.

7. Once you are happy:

    - Publish the release tarball (currently, we use box.com Berkeley account).
    - Update the Download section of the BerkeleyGW website.
    - Commit (locally) any changes you may have made.
    - Tag the current version, e.g.: `git tag -a v2.1 -m 'version 2.1'`.
       You may need to add a `--force` to the command if you want to correct
       a wrong tag.
    - Push these changes directly to the `upstream` and `origin` repos:
      `git push --tags -f origin master && git push --tags -f upstream master`
    - Upload a new version of the manual. Read the instructions how to create
      and publish the user manual (`documentation/user_manual/README.md`) from
      within the development version of BerkeleyGW. However, perform the actual
      steps to assemble the manual from within the released version of
      BerkeleyGW, to avoid publishing development tags/features.
    - If necessary, change the default page for the manual. For instance, if you
      released version 2.1 of the code, you want <http://manual.berkeleygw.org> to
      redirect to <http://manual.berkeleygw.org/2.1>. To do that:
        - Go to the AWS console
        - Click on Services, then S3
        - Click on `manual.berkeley.org` (on the name, not the check box)
        - Click on `index.html`
        - Click on the tab `Properties`, then click `Metadata`
        - Click on the box next to `Website-Redirect-Location`, click on `Edit`,
          and put the desired version to redirect to, e.g., `/2.1/`.
        - Save.
