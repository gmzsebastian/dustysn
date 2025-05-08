.. _usage:

Usage
=====

Basic Usage
-----------

The simplest way to use ``dustysn`` is through the main function:

.. code-block:: python

    import dustysn

    # Fit dust

This function will produce a plot with all diagnostics to determine if the transient is nuclear.

Requirements
------------

For ``dustysn`` to work properly:

* The SED must be...

Core Functions
--------------

fit_data
~~~~~~~~
The ``fit_data`` function fits the data:

.. code-block:: python

    from dustysn.utils import fit_data
    
    fit_data("phtometry")

