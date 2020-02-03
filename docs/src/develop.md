# Developers Guide
`InfiniteOpt` is large project with a lot of opportunity for development. As
such we warmly welcome and encourage contributions. This page serves as the guide
of how contributions can be made and how we prefer that they be carried out.

## Contribution Roadmap
This section will provide a birds-eye view on how to make a contribution to this
project. More precise details such as the preferred style are detailed in the
sections further below.

So you want to help improve `InfiniteOpt`, awesome and thank you! Let's walk
through step by step how this should be done.

  1. Setup a GitHub account if you do not already have one. Here is the
     [link](https://github.com/join) to do so.
  2. Select a specific task to develop that is well-defined. This can as simple as
     correcting/clarifying a documentation page or as involved as implementing a
     more efficient data management paradigm. With a task in mind, please start a new
     issue [here](https://github.com/pulsipher/InfiniteOpt.jl/issues) in the
     `InfiniteOpt` repository. Also, this a good place find tasks to contribute to by
     browsing what open issues are (especially ones with the tag `good first issue`).
     Note that if your proposed contribution corresponds to an existing issue please
     do not make a new issue. A guide to using issues in GitHub is located
     [here](https://guides.github.com/features/issues/).
  3. Fork the `InfiniteOpt` repository to your GitHuB account. Only core
     developers have permissions to modify `InfiniteOpt` directly, thus other need
     to fork it which essentially amounts to creating their own linked copy. This is
     done by clicking the `Fork` button at the top left corner on the main repository
     page [here](https://github.com/pulsipher/InfiniteOpt.jl).
  4. Install Git on your computer. Git is an open source version control program
     for repositories (it is why GitHub uses the word Git). This needed to manipulate
     the repository (all the package files) locally on your computer. A simple google
     search should indicate how his should be done for your computer. I personally
     prefer [Git for Windows](https://gitforwindows.org/) as a Windows user.
  5. Now you need to install your forked version of `InfiniteOpt` in Julia on your
     computer. This needs to be done via the `dev` command in the package manager
     so you can edit it. The syntax is as follows:
     ```julia
     (v1.3) pkg> dev https://github.com/username-here/InfiniteOpt.jl
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
  8. Create a pull request. Go [here](https://github.com/pulsipher/InfiniteOpt.jl)
     to `InfiniteOpt`'s main page and create a pull request drawing from your forked
     repository. A step by step explanation is provided
     [here](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/creating-a-pull-request-from-a-fork).
  9. Make necessary changes if the tests fail and/or we ask you to make specific
     changes. The Codecov tests will ensure every new line of code is tested at least
     once with the new test functions and Travis CI will ensure that the tests pass
     on a range of operating systems and Julia versions.
  10. That's it. Once the new additions are ready, we will merge them into the
      main repository.
  11. Contribute more by repeating steps 2 and 6-10. Just make sure to update your
     forked repository before getting started which can be done as explained
     [here](https://github.com/KirstieJane/STEMMRoleModels/wiki/Syncing-your-fork-to-the-original-repository-via-the-browser).
     Also, be sure to pull the updated repository unto your computer before getting
     started.

## Style Guide
Below we detail the formatting, naming, and organizational styles used in
`InfiniteOpt`. We kindly ask developers to adhere to these practices in efforts
to foster consistency.

### File Organization
Files for `InfiniteOpt` are principally stored in 5 locations:
 - the base directory `./`,
 - the source file directory `./src/`,
 - the source code testing directory `./test/`,
 - the documentation source file directory `./docs/`,
 - and the example scripts directory `./examples/`.

The base directory is for certain critical package files such as the `README.md`
file and the CI (virtual testing service) configuration files. Files should NOT
be added or removed from this directory, but may be modified as needed.

Naturally, the source directory is where all the source code files are located.
The principal file `InfiniteOpt.jl` is where the main module is defined, all
source code files are included, and all methods/datatypes/macros are exported. This
file shouldn't contain any function or datatype definitions directly, but rather
should include source files containing such via `include("file_name.jl")`. Where
possible new datatypes should be defined in `datatypes.jl` and new methods should
be defined in the appropriate file (e.g., a new parameter method should be
defined in `parameters.jl`). New files can be added as necessary to help with
organization and to prevent a particular file from becoming too long. Also, note
that any submodule (e.g., `TranscriptionOpt`) should be defined within its own
sub-directory named after itself.

The test directory contains all the files in appropriate organization to test
all of the methods, datatypes, and macros defined in the source files. The file
structure here should emulate that of the `./src/` directory since each file
should by systematically tested as described be below in the Unit Tests section.
Here the principle file is `runtests.jl` which serves as the backbone for all the
unit testing. Again, no explicit tests should be contain in it, but rather
inclusions of test files via `include("file_name")`.

The documentation directory follows a particular structure as explained in the
documentation for `Documenter.jl`. Here the root directory `./docs/` contains
`make.jl` which is the script that generates the documentation via `Documenter.jl`.
The `Project.toml` includes the packages necessary to do this. The `./docs/src/`
sub-directory is where source code is stored to build the documentation pages.
When building the documentation locally, a `./docs/build` directory will also
appear that stores the built HTML files. However, this directory is not tracked
by Git and any changes here will be ignored.

The example directory contains scripted use examples of `InfiniteOpt`. Each
example should be stored in single `.jl` file where possible. However, other more
complex examples that use multiple files should be stored in an appropriately
named folder.

Please note that all file/folder names should use complete names and avoid
abbreviations where possible unless the abbreviations are unambiguous and common
knowledge. Moreover, names should be lowercase and use underscores between
words: `example_file_name.jl`.  

### Julia Code
Here we detail the programmatic style used with `InfiniteOpt`. This is done in an
effort to make this package intuitive for new-comers and to ease development. This
style closely follows that of `JuMP.jl` with similar deviations from typical Julia
styles. Please refer to the `JuMP` style guide
[here](http://www.juliaopt.org/JuMP.jl/stable/style/) as this is the style used
by `InfiniteOpt`.

In addition, we adopt the following practices:
 - All names should be meaningful and readily identifiable. This is bad:
   ```julia
   x = y2 - cp
   ```
   This is good:
   ```julia
   new_pizza_cost = old_pizza_cost - discount
   ```
   This will make lines longer, but much more understandable.
 - Avoid the use explicit numeric values (i.e., magic numbers):
   This is bad:
   ```julia
   tax = 0.07 * total_price
   ```
   Typically, this will employ the use of constants via `const`
   This is good:
   ```julia
   const TAX_RATE = 0.07
   tax = TAX_RATE * total_price
   ```
   Exceptions to this rule include the use of `1`, `1.0`, `0`, `0.0`, `-1`,
   `-1.0`, `Inf`, and `-Inf`.
 - Where possible use `eachindex` to iterate over an datatype:
   This is bad:
   ```julia
   for i in 1:length(A)
       A[i] = i
   end
   ```
   This is good:
   ```julia
   for i in eachindex(A)
       A[i] = i
   end
   ```
 - All function arguments and struct elements should be typed. Also, function
   outputs should be typed where possible unless nothing is returned.   
   This is bad:
   ```julia
   function my_new_function(arg1, arg2)
       return arg1 + arg2
   end
   struct MyNewStruct
       thing1
       thing2
   end
   ```
   This is good:
   ```julia
   function my_new_function(arg1::Int, arg2::Int)::Int
       return arg1 + arg2
   end
   struct MyNewStruct
       thing1::Int
       thing2::String
   end
   ```
 - Type dispatch should be used instead of conditional statements based on type:
   This is bad:
   ```julia
   function my_new_function(arg::AbstractType)
     if arg isa Type1
         temp = arg + 1
     elseif arg isa Type2
         temp = 0
     end
     # do more stuff with temp
     return temp
   end
   ```
   This is good:
   ```julia
   ## Internal dispatch for my_new_function
   # Type1
   function _my_internal_function(arg::Type1)
       return arg + 1
   end
   # Type2
   function _my_internal_function(arg::Type1)
       return 0
   end
   # Fallback
   function _my_internal_function(arg::AbstractType)
       error("Unrecognized type...")
       return
   end
   # Main method
   function my_new_function(arg::AbstractType)
     temp = _my_internal_function(arg)
     # do more stuff with temp
     return temp
   end
   ```
 - Functions should be built in a modular manner to avoid code repetition and
   excessively long function definitions.

In addition to the above guidelines, contributions should be structured such that
extensions are readily possible without having to rewrite all of the associated
functions. The ability to easily facilitate extensions is a core goal of
`InfiniteOpt` and this should be kept in mind when developing contributions.

TODO add example.

### Docstrings and Comments
Here we discuss the use of Docstrings and comments in `InfiniteOpt`. All public
functions, macros, and datatypes should have a Docstring. This is enables the
help query tool in Julia and is needed for inclusion in the documentation pages.
For functions and macros the format should follow the form:
````julia
"""
    my_new_function(arg1::Type, [arg2::Type = 0; karg1::Type = 42])::Type

Precise and concise description of what `my_new_function` does and what it
returns (also what will cause it will trigger errors). This is in markdown
format.

**Example**
```jldoctest; setup = :(using JuMP, InfiniteOpt; [other commands])
julia> my_new_function(input...)
expected_output
```
"""
function my_new_function(arg1::Type, arg2::Type = 0; karg1::Type = 42)::Type
    return arg1 + arg2 + karg1
end
````
Notice that the function is declared at the top with an ident and the optional
arguments are enclosed within square brackets. This can be spaced over several
lines if there are too many arguments to fit on one line. Also, the `jldoctest`
should be used where possible (this actually runs and tests that the example
works), but when this is not possible it can be replaced with `julia` which
will provide syntax highlighting without running the code.

For datatypes Docstrings should follow the form:
```julia
"""
    MyNewStruct

Precise and concise description of what this is.

**Fields**
- `element1::Type` Description of what this is.
- `element2::Type` Description of what this is.
"""
struct MyNewStruct
    element1::Type
    element2::Type
end
```
Note that if the struct is parametric and/or has inheritance, this information
should also be shown in the header. For example, we have that
`InfOptParameter{T <: AbstractInfiniteSet} <: JuMP.AbstractVariable`.

For more docstring information please visit the Julia documentation
[here](https://docs.julialang.org/en/v1/manual/documentation/index.html).

Furthermore, all internal functions and datatypes should have an appropriate
commented description of what they do above them. This should follow the
form:
```julia
# Description of what _my_internal_function does. Bla Bla Bla Bla Bla Bla Bla
# Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla Bla.
function _my_internal_function(arg1::Type, arg2::Type)::Type
    return arg1 + arg2
end
```

Finally, we encourage a healthy usage of comments throughout source code to
enhance its readability. A simple comment before a complex block of code can make
all the difference.

### Unit Tests
A nice attribute of `InfiniteOpt` is that it is near perfect code testing
coverage. This success is due to strictly testing every method and macro rigorously
such that every line is called. This has been very advantageous in detecting
many bugs which can be difficult to anticipate given the quantity of source code.
Thus, tests must be created/updated to cover any new additions/changes in the
`./src/` directory.

The `runtests.jl` file serves as the principal backbone for doing this. We use
a nested `@testset` structure using `Test.jl`. Please refer to the documentation
[here](https://docs.julialang.org/en/v1/stdlib/Test/) to learn about the relevant
testing macros. The structure typically groups related functions together where
each function/macro/datatype is tested via a `@testset` that employs a number of
tests via `@test` that thoroughly test it. This is typically of the form:
```julia
@testset "my_new_function" begin
    @test my_new_function(input1) == expected_output1
    @test my_new_function(input2) == expected_output2
    @test my_new_function(input3) == expected_output3
    .
    .
    .
end
```

Thus, a function's `@testset` should be updated when the respective function has
been modified. Moreover, a new `@testset` should be added for each new
function/macro. New function tests should be implemented in an order such that
any other functions/macros they depend on are tested first.

Also, where possible please include comments to explain what is going on.

Please refer to `InfiniteOpt/test/` for examples.

### Documentation Pages
Documentation in `InfiniteOpt` is generated via
[`Documenter.jl`](https://github.com/JuliaDocs/Documenter.jl). Please refer to
its documentation to learn about how to use it.

The source markdown files stored in `./docs/src/` are what comprise the source
code for the documentation pages and are principally what should be updated. A
guide for markdown syntax is provided [here](https://www.markdownguide.org/).
Also, note that `Documenter` enables unique functionality in addition to this
general guide.

When a new Docstring is created as described above, it should be included on the
appropriate guide page in the `@docs` block at the bottom. Moreover, content
should be added in an appropriate section above (or perhaps in a new section)
that overviews how to implement the new capabilities in an example driven fashion.
These examples should use `jldoctest`s where possible as well to assess whether
the example code is functional.

Documentation content should be concise and use examples and lists where possible
to provide a more visual guide. Also, we ask that passive voice be avoided.

Be sure to test the documentation first locally by running `make.jl` to check
for problems which may include:
 - unrecognized docstrings
 - failed doctests
 - faulty links
 - unrecognized formats
 - missing package dependencies
 - etc.

!!! note
    Doctests should be setup to pass on Linux and not other operating systems
    since Linux is used to build the website. For example, in the symbol `∀` is
    used in Linux while in Windows the phrase `for all` is used. Thus, a number
    of doctests will fail when run on Windows which is fine. However, users
    should choose the expected output to be what Linux will return.
