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

During the re-base, you may see some conflicts and you need to resolve them. You will need to resolve the conflict one-by-one. For each conflict, you can open VS Code and click "Source Control" in the left panel. Then, you can see the conflicting files with an escalation mark. Open these files, and you will see the code with the conflict marked in colors. You will most likely choose "Accept Current Change", which will use the latest code from v2512-orig (you might need to slightly tweak the code to fix any potential AD compilation errors). Once you resolve the conflicts, run `git add .` to add the edited files, and then run `git rebase --continue` to resolve the next conflict. Right after you run `git rebase --continue`, you will see some pop-up texts about the conflict, and you need to close the window by running `:wq` (assuming you use vim). You need to repeat this for all conflicts!

After ALL the conflicts are solved, you can do

<pre>
git push -u origin v2512-ad
</pre>

