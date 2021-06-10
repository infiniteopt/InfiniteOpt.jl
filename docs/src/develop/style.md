# Style Guide
Below we detail the formatting, naming, and organizational styles used in 
`InfiniteOpt`. We kindly ask developers to adhere to these practices in efforts 
to foster consistency.

## File Organization
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

## Julia Code
Here we detail the programmatic style used with `InfiniteOpt`. This is done in an 
effort to make this package intuitive for new-comers and to ease development. This 
style closely follows that of `JuMP.jl` with similar deviations from typical Julia 
styles. Please refer to the  
[`JuMP` style guide](https://jump.dev/JuMP.jl/v0.21.8/developers/style/) as this 
is the style used by `InfiniteOpt`.

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
 - All function arguments and `struct` elements should be typed. Also, function
   outputs should be typed where possible.   
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
   function my_new_function(arg::AbstractType)::ReturnType
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
   function _my_internal_function(arg::Type1)::Int
       return arg + 1
   end
   # Type2
   function _my_internal_function(arg::Type1)::Int
       return 0
   end
   # Fallback
   function _my_internal_function(arg::AbstractType)
       error("Unrecognized type...")
   end
   # Main method
   function my_new_function(arg::AbstractType)::ReturnType
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

## Docstrings and Comments
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
```julia-repl
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
lines if there are too many arguments to fit on one line. 

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

## Unit Tests
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

## Documentation Pages
Documentation in `InfiniteOpt` is generated via 
[`Documenter.jl`](https://github.com/JuliaDocs/Documenter.jl). Please refer to 
its documentation to learn about how to use it.

The source markdown files stored in `./docs/src/` are what comprise the source 
code for the documentation pages and are principally what should be updated. A 
guide for markdown syntax is provided [here](https://www.markdownguide.org/). 
Also, note that `Documenter` enables unique functionality in addition to this 
general guide.

When a new Docstring is created as described above, it should be included in the 
appropriate `@docs` block on its corresponding manual page. Moreover, content 
should be added in an appropriate section (or perhaps in a new section) in the 
guide that overviews how to implement the new capabilities in an example driven 
fashion. These examples should use `jldoctest`s where possible as well to assess 
whether the example code is functional.

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

## Case Study Examples
We use [`Literate.jl`](https://fredrikekre.github.io/Literate.jl/v2.8/) to run 
the case studies in `./docs/src/examples/` and generate markdown files that are 
incorporated into the documentation for the `Examples` sections. 

A new case study example can be added to an appropriate sub-folder of 
`./docs/src/examples/` (a new sub-folder can be made if needed). The example file 
should be `.jl` file that uses comments in accordance with `Literate.jl`'s format. 
This is exemplified below:
```julia
# # My Example Name
# Text to introduce my example...

# ## Background 
# Text that describes the problem we are trying to solve. We can also include 
# latex math like ``x^2`` and math blocks such as:
# ```math 
# x^2 + y = 1
# ```

# ## Formulation
# Text to introduce as needed...
using InfiniteOpt, Clp # import the needed packages

## This comment type will be part of the code block
model = InfiniteModel(Clp.Optimizer) # add side comments to code

# This comment type will be Markdown again, thus breaking up the code block
@infinite_parameter(model, t in [0, 1], num_supports = 42)
@variable(model, y >= 0, Infinite(t))

optimize!(model)

## TODO add more code

# ### Maintenance Tests
# These are here to ensure this example stays up to date. 
using Test
@test termination_status(model) == MOI.OPTIMAL
@test has_values(model)
## Add more tests as appropriate
```
The above file will then be tested and incorporated into the documentation when 
`./docs/make.jl` is called. Notice that we also have the "Maintenance Test" 
section at the end that will be used to run checks to ensure the example script 
is working as expected (this helps ensure the documentation is up to date).

The above example would produce the following markdown file via `Literate.jl`:
````markdown
# My Example Name
Text to introduce my example...

## Background
Text that describes the problem we are trying to solve. We can also include
latex math like ``x^2`` and math blocks such as:
```math
x^2 + y = 1
```

## Formulation
Text to introduce as needed...

```julia
using InfiniteOpt, Clp # import the needed packages

# This comment type will be part of the code block
model = InfiniteModel(Clp.Optimizer) # add side comments to code

```

```julia
An InfiniteOpt Model
Feasibility problem with:
Finite Parameters: 0
Infinite Parameters: 0
Variables: 0
Derivatives: 0
Measures: 0
Optimizer model backend information:
Model mode: AUTOMATIC
CachingOptimizer state: EMPTY_OPTIMIZER
Solver name: Clp
```

This comment type will be Markdown again, thus breaking up the code block

```julia
@infinite_parameter(model, t in [0, 1], num_supports = 42)
@variable(model, y >= 0, Infinite(t))

optimize!(model)

# TODO add more code
```

### Maintenance Tests
These are here to ensure this example stays up to date.

```julia
using Test
@test termination_status(model) == MOI.OPTIMAL
@test has_values(model)
# Add more tests as appropriate
```

```
Test Passed
```

---

*This page was generated using [Literate.jl](https://github.com/fredrikekre/Literate.jl).*


````
