# Visual QA contract

Formal figures pass both programmatic and visual checks.

Programmatic checks: all required files exist; PNG DPI and dimensions match the
contract; PDF/SVG exports exist or have a sidecar reason; fonts and color maps
are recorded; source CSV and input hashes exist; no NaN/Inf is plotted without
an explicit missing label; panel order and sort keys are recorded.

Visual checks at final size: no clipping or overlap; all tick labels and legends
are readable; panel labels align; dense points do not hide the intended trend;
gray/empty `NOT_ESTIMABLE` slots remain visible; categorical colors are stable;
the figure remains interpretable in grayscale; the caption does not exceed the
evidence ceiling. Re-render after every correction.

This QA does not validate the scientific model. It validates that the rendered
artifact faithfully exposes the declared model and its limitations.

