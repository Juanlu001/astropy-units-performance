## Setup

* astropy
* poliastro (optional, for some examples)
* ipython (recommended)
* line_profiler

## Notes

Some critical functions:

* Quantity.is_unit()
* Quantity.__init__
* CompositeUnit.__init__

Ideas:

* Caching common composite units would help reducing number of operations
* Array * Unit is still slow
