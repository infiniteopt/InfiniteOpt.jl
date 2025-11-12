# [Extensions](@id ext_manual)
A technical manual for extension packages natively hosted by `InfiniteOpt`. See the 
respective [guide](@ref ext_docs) for more information.

## InfiniteInterpolations
Enabled via `import InfiniteOpt, Interpolations`.
```@docs
JuMP.value(::Union{GeneralVariableRef,JuMP.GenericAffExpr{Float64, GeneralVariableRef},JuMP.GenericQuadExpr{Float64, GeneralVariableRef},JuMP.GenericNonlinearExpr{GeneralVariableRef},InfOptConstraintRef}, ::Interpolations.InterpolationType)
```

## InfiniteMathOptAI
Enabled via `import InfiniteOpt, MathOptAI`.
```@docs
MathOptAI.add_variables
```