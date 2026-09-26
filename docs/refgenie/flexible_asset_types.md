# Refgenie Recipe and Asset Class System

## Overview

Refgenie’s new recipe and asset class system introduces a flexible, extensible, and user-driven approach to defining and managing reference genome assets. This system eliminates the previous reliance on internally-defined asset types, empowering users and tool developers to define, share, and distribute their own asset types and recipes.

## Key Concepts

### Asset Classes

- **Asset classes** define the structure and seek keys (files or directories) that make up an asset type (e.g., a FASTA, GTF, or index).
- Asset classes are no longer hardcoded in refgenie. Instead, they are defined in external YAML or JSON files, which can be created, modified, and shared by anyone.
- This means new asset types can be introduced without modifying refgenie’s source code.

### Recipes

- **Recipes** describe how to build an asset of a given class from input assets.
- Recipes are also defined externally and can be distributed independently of refgenie itself.
- Recipes specify the required input asset classes, parameters, and the steps to build the output asset.

### Data Channels

- Both asset classes and recipes can be distributed through **data channels**—remote or local repositories that refgenie can subscribe to.
- This enables community-driven sharing and rapid adoption of new asset types and build methods.

## Benefits

### Decoupling and Extensibility

- **No more internal asset type lock-in:** Refgenie no longer restricts users to a fixed set of asset types. Anyone can define new asset classes and recipes.
- **Community-driven innovation:** Tool developers and users can publish and share new asset types and build recipes, fostering a collaborative ecosystem.
- **Distribution via data channels:** Asset classes and recipes can be versioned and distributed through data channels, making it easy to adopt new standards or methods.

### Solving the 1:1 Asset-Recipe Coupling Problem

Previously, recipes were tightly coupled to asset types. For example, a genome could only have one `fasta` asset (e.g., `hg38/fasta`), and all recipes requiring a FASTA as input would use this single asset. This design made it impossible to have multiple assets of the same type under a genome namespace (e.g., both `hg38/fasta` and `hg38/fasta_txome`).

**With the new system:**

- Multiple assets of the same class can coexist under a genome (e.g., `hg38/fasta`, `hg38/fasta_txome`, etc.).
- Recipes can specify which asset of a given class to use as input, allowing for more complex and flexible workflows.
- This decoupling enables scenarios where, for example, you can build a transcriptome FASTA (`fasta_txome`) alongside the primary genome FASTA, and use either as input to downstream recipes.

## Example Workflow

1. **Define an asset class:**  
   Create a YAML file describing the structure of a new asset type (e.g., a custom index).
2. **Write a recipe:**  
   Create a recipe file specifying how to build this asset from input assets.
3. **Distribute via data channel:**  
   Publish the asset class and recipe to a data channel for others to use.
4. **Build assets:**  
   Users can now build multiple assets of the same class under a genome, and recipes can consume any of these as inputs.

## Conclusion

The recipe and asset class system in refgenie makes the platform more modular, extensible, and community-driven. By decoupling asset types from internal definitions and enabling external, user-defined recipes and asset classes, refgenie supports a much broader range of use cases and workflows.
