How to AD OpennFOAM
===================

First, download this repository and add OpenFOAM's original GitLab repo as the upstream (NOTE: this needs to be done only once):

<pre>
git clone https://github.com/DAFoam/OpenFOAM-AD.git && \
cd OpenFOAM-AD && \
git remote add upstream https://gitlab.com/openfoam/core/openfoam.git && \
git push -u origin --all && \
git push origin --tags
</pre>

Then, copy OpenFOAM's existing tags to create the orig branch as the base. Here we create a `v2412-orig` branch based on OpenFOAM's `OpenFOAM-v2412` tag. Don't change the v2412-orig branch as it is used as a reference.

<pre>
git fetch upstream --tags && \
git checkout -b v2412-orig OpenFOAM-v2412 && \
git push -u origin v2412-orig
</pre>

Next, we can create an ad branch based on the orig branch.

<pre>
git checkout -b v2412-ad v2412-orig
</pre>

Now, we can add and push our AD implementations to the v2412-ad branch.

When there is a new tag from OpenFOAM's GitLab repo, e.g., OpenFOAM-v2512, we can re-base the v2412-ad to it by first creating a reference v2512 branch called v2512-orig

<pre>
git fetch upstream --tags && \
git checkout -b v2512-orig OpenFOAM-v2512 && \
git push -u origin v2512-orig
</pre>

Then, we can do the re-base

<pre>
git checkout v2412-ad && \
git branch v2512-ad && \
git checkout v2512-ad && \
git rebase v2512-orig
</pre>

During the re-base, you may see some conflicts and you need to resolve them before you can push it to the new v2512-ad branch. After the conflicts are solved, you can do

<pre>
git push -u origin v2512-ad
</pre>

