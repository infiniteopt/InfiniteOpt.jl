"""
    add_measure(model::InfiniteModel, v::Measure)
Add a measure to the `InfiniteModel` object in an analagous way to `JuMP.add_variable`.
"""
function add_measure(m::InfiniteModel, meas::Measure)
    m.next_meas_index += 1
    mref = MeasureRef(m, m.next_meas_index)
    m.measures[mref.index] = meas
    mref.model.meas_to_name[mref.index] = _make_meas_name(meas)
    return mref
end

# Set a default weight function
_w(t) = 1

# Define constructor functions
DiscretizeData() = DiscretizeData(_w, [], [])
DiscretizeData(coeffs::Vector, pts::Array) = DiscretizeData(_w, coeffs, pts)
Measure(func::InfiniteExpr, data::AbstractMeasureData) = Measure("measure", func, data)

# Method for defining method on the fly
function measure(name::String, expr::InfiniteExpr, data::AbstractMeasureData)
    meas = Measure(name, expr, data)
    model = _get_model_from_expr(expr)
    return add_measure(model, meas)
end

# Method for defining method on the fly
function measure(expr::InfiniteExpr, data::AbstractMeasureData)
    meas = Measure(expr, data)
    model = _get_model_from_expr(expr)
    return add_measure(model, meas)
end

# Parse the model pertaining to an expression
function _get_model_from_expr(expr::JuMP.AbstractJuMPScalar)
    if expr isa JuMP.AbstractVariableRef
        return expr.model
    elseif expr isa JuMP.GenericAffExpr
        aff_vars = [k for k in keys(expr.terms)]
        if length(aff_vars) > 0
            return aff_vars[1].model
        else
            return
        end
    elseif expr isa JuMP.GenericQuadExpr
        aff_vars = [k for k in keys(expr.aff.terms)]
        if length(aff_vars) > 0
            return aff_vars[1].model
        else
            var_pairs = [k for k in keys(expr.terms)]
            if length(var_pairs) > 0
                return var_pairs[1].a.model
            else
                return
            end
        end
    else
        return expr.m
    end
end

# Parse the string for displaying a measure
function _make_meas_name(meas::Measure)
    return string(meas.name, "(", JuMP.function_string(JuMP.REPLMode, meas.func), ")")
end

"""
    JuMP.name(mref::MeasureRef)
Extend the `JuMP.name` function to accomodate measure references.
"""
JuMP.name(mref::MeasureRef) = mref.model.meas_to_name[mref.index]

# TODO Add manipulation functions like variables have
