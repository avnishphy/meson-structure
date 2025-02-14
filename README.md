# meson-structure
Meson structure analysis

```mermaid
flowchart TD
    A[HCAL Clusters from ZDC]
    B[FarForwardNeutralsReconstruction]
    C{Cluster Classification}
    D[Gamma-like Clusters]
    E[Non-Gamma Clusters]
    F[Reconstructed Gamma Candidates]
    G[Accumulate Energy for Neutron Candidate]
    H[Reconstructed Neutron Candidate]
    I[Inputs: 1 Neutron & 2 Gammas]
    J[FarForwardLambdaReconstruction]
    K[Iterative Vertex Finding]
    L[Lambda Candidate (n + π⁰)]
    M[Decay Products: n, γ, γ]
    N[Quality Check: Mass within Tolerance]
    O[Final Lambda & Decay Product Collection]

    A --> B
    B --> C
    C -- Pass Cuts --> D
    C -- Fail Cuts --> E
    D --> F
    E --> G
    G --> H
    F & H --> I
    I --> J
    J --> K
    K --> L
    L --> N
    N -- Pass --> M
    M --> O
```

```mermaid
flowchart TD
    A[FarForwardNeutralsReconstruction_factory]
    B[FarForwardNeutralsReconstruction Module]
    C[Processes HCAL Clusters]
    D[Classifies Clusters into Gammas & Neutron]
    E[Applies Energy Corrections]
    
    F[FarForwardLambdaReconstruction_factory]
    G[FarForwardLambdaReconstruction Module]
    H[Receives Neutron & Gamma Candidates]
    I[Applies Coordinate Rotation]
    J[Performs Iterative Vertex Finding]
    K[Validates Lambda Invariant Mass]
    L[Rotates Back & Boosts Decay Products]
    
    M[Global EICRecon Pipeline (reco.cc)]
    
    A --> B
    B --> C
    C --> D
    D --> E
    E --> M
    
    F --> G
    G --> H
    H --> I
    I --> J
    J --> K
    K --> L
    L --> M

```