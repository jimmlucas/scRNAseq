# Normalization Approval Phase Test

This test validates the `VALIDATE_NORMALIZATION_APPROVAL` control gate.

It checks that a reviewed normalization policy can only pass when human
approval metadata, scientific identity, strategy parameters and application
state are valid. The phase does not run normalization and does not modify the
Seurat object.

