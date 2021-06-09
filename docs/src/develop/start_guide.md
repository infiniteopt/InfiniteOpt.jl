# [Getting Started](@id contribute_guide)
`InfiniteOpt` is a large project with a lot of opportunity for development. As 
such we warmly welcome and encourage contributions. This page serves as the guide 
of how to start contributing.

Before starting please review our 
[Code of Conduct](https://github.com/pulsipher/InfiniteOpt.jl/blob/master/CODE_OF_CONDUCT.md).

## Step-by-Step
This section will provide a birds-eye view on how to make a contribution to this 
project.

So you want to help improve `InfiniteOpt`, awesome and thank you! Let's walk 
through step by step how this should be done.

  1. Setup a GitHub account if you do not already have one. Here is the 
     [link](https://github.com/join) to do so.
  2. Select a specific task to develop that is well-defined. This can as simple as 
     correcting/clarifying a documentation page or as involved as implementing a 
     more efficient data management paradigm. With a task in mind, please start a new 
     issue [here](https://github.com/pulsipher/InfiniteOpt.jl/issues) in the 
     `InfiniteOpt` repository. Also, this is a good place to find tasks to contribute to by 
     browsing what open issues are (especially ones with the tag `good first issue`). 
     Note that if your proposed contribution corresponds to an existing issue please 
     do not make a new issue. A guide to using issues in GitHub is located 
     [here](https://guides.github.com/features/issues/).
  3. Fork the `InfiniteOpt` repository to your GitHub account. Only core 
     developers have permissions to modify `InfiniteOpt` directly, thus others need 
     to fork it which essentially amounts to creating their own linked copy. This is 
     done by clicking the `Fork` button at the top left corner on the main repository 
     page [here](https://github.com/pulsipher/InfiniteOpt.jl).
  4. Install Git on your computer. Git is an open source version control program 
     for repositories (it is why GitHub uses the word Git). This is needed to manipulate 
     the repository (all the package files) locally on your computer. A simple Google 
     search should indicate how this should be done for your computer. I personally 
     prefer [Git for Windows](https://gitforwindows.org/) as a Windows user. 
  5. Now you need to install your forked version of `InfiniteOpt` in Julia on your 
     computer. This needs to be done via the `dev` command in the package manager 
     so you can edit it. The syntax is as follows:
     ```julia
     (v1.6) pkg> dev https://github.com/username-here/InfiniteOpt.jl
     ```
     We also recommend you install [`Revise.jl`](https://github.com/timholy/Revise.jl) 
     which is very useful when developing packages in Julia.
  6. Develop your contribution. Please follow the style guides featured in the 
     sections below. A programmatic contribution will involve the following parts:
      - editing/adding code to the `.jl` files in the `src` (source) directory
      - adding a docstring for each public function/datatype
      - including comments that describe each internal function/datatype
      - adding a unit-testing for each function in the appropriate test files in the 
        `test` directory
      - adding documentation of the new functionality in the appropriate place in the 
        documentation by adding the files in the `docs/src` directory.
     These aspects are detailed further in the sections below.
  7. Commit and push your changes to your forked repository. This is done via Git 
     using your preferred interface and one should pull, add, commit, and then push 
     the changes. Using a bash terminal it would look like this:
     ```bash
     username@ubuntu:~$ cd repo_directory
     username@ubuntu:~/repo_directory$ git pull origin master
     username@ubuntu:~/repo_directory$ git add *
     username@ubuntu:~/repo_directory$ git commit -m "insert commit message here"
     username@ubuntu:~/repo_directory$ git push origin master
     ```
     We recommend using [VsCode](https://www.julia-vscode.org/) as an editor and 
     as a GUI for interfacing with Git.
  8. Create a pull request. Go [here](https://github.com/pulsipher/InfiniteOpt.jl) 
     to `InfiniteOpt`'s main page and create a pull request drawing from your forked 
     repository. A step by step explanation is provided 
     [here](https://docs.github.com/en/github/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request-from-a-fork).
  9. Make necessary changes if the tests fail and/or we ask you to make specific 
     changes. The Codecov tests will ensure every new line of code is tested at least 
     once with the new test functions and the GitHub Actions CI will ensure that 
     the tests pass on a range of operating systems and Julia versions.
  10. That's it. Once the new additions are ready, we will merge them into the 
      main repository.
  11. Contribute more by repeating steps 2 and 6-10. Just make sure to update your 
     forked repository before getting started which can be done as explained 
     [here](https://github.com/KirstieJane/STEMMRoleModels/wiki/Syncing-your-fork-to-the-original-repository-via-the-browser). 
     Also, be sure to pull the updated repository unto your computer before getting 
     started.
