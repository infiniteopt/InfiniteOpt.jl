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
MeasureData() = MeasureData(_w, [], [])
MeasureData(coeffs::Vector, pts::Array) = MeasureData(_w, coeffs, pts)
Measure(func::InfiniteExpr, data::MeasureData) = Measure("measure", func, data)

# Method for defining method on the fly
function measure(name::String, expr::InfiniteExpr, data::MeasureData)
    meas = Measure(name, expr, data)
    if expr isa InfiniteVariableRef
        model = expr.model
    elseif expr isa InfiniteAffExpr
        model = [k for k in keys(expr.terms)][1].model
    else
        aff_vars = [k for k in keys(expr.terms)]
        if length(aff_vars) > 0
            model = aff_vars[1].model
        else
            model = [k for k in keys(expr.terms)][1].a.model
        end
    end
    return add_measure(model, meas)
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
