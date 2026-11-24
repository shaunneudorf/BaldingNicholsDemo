# Balding-Nichols likelihood ratio demo

A Blazor WebAssembly page that shows how the Balding-Nichols sampling formula
builds a likelihood ratio for a suspect, with or without an assumed contributor
(such as a victim) in a two-person mixture.

The sampling formula used for every draw is

    P(next allele = A) = [x·θ + (1−θ)·Pr(A)] / [1 + (n−1)·θ]

where x is the number of copies of A already seen and n is the number of alleles
already seen. Deleting the victim reduces the result to the NRC II 4.2 match
probabilities, and θ = 0 gives the product rule. Unit tests cover both.

## Allele frequencies

`src/BaldingNichols.Web/wwwroot/data/Caucasian_database.csv` holds the frequency tables. Pick a locus and the
suspect's allele frequencies fill in automatically; any value can be edited, and "Use database value" restores it.

- First column is the allele designation, the other columns are loci.
- A row labelled `N` gives the number of alleles in the database for each locus.
- Balding-Nichols draws use `(x + 1/k) / (N + 1)` for every allele, where x is the number of copies in the
  database, k is the number of designations with non-zero observations and N is the number of alleles. For an allele
  not observed (x = 0) this is `(1/k) / (N + 1)`.
- The product rule uses the database frequency directly, and `5/2N` for an allele not observed. In this file N is
  already the allele count (that is, 2N for a count of individuals), so 5/2N is computed as 5 / N. The 5 is the
  constant `ProductRuleUnseenNumerator` in `AlleleDatabase.cs`.
- The file holds frequencies rather than counts, so x is recovered as frequency × N and rounded.

To use a different database, replace the CSV with the same layout. Changing the file name also means changing
the path in `Calculator.razor` and the link in the test project.

## Layout

- `src/BaldingNichols.Core` – calculator and frequency database reader (no UI dependencies)
- `src/BaldingNichols.Web` – Blazor WebAssembly front end
- `tests/BaldingNichols.Tests` – xUnit tests

## Run locally

Open `BaldingNicholsDemo.sln` in Visual Studio (2022 17.8+ with .NET 8), set
`BaldingNichols.Web` as the startup project and press F5. Or:

    dotnet run --project src/BaldingNichols.Web

Run the tests with `dotnet test`.

## Host on GitHub Pages

1. Push this folder to a GitHub repository whose default branch is `main`.
2. In the repository, go to Settings > Pages and set Source to **GitHub Actions**.
3. Push to `main`. The workflow in `.github/workflows/deploy.yml` runs the tests,
   publishes the app and deploys it to `https://<user>.github.io/<repo>/`.

If the repository is named `<user>.github.io`, the site lives at the root, so change the
`sed` line in the workflow to leave `<base href="/" />` unchanged.

## Limits of the model

Pr(E|Hp) is taken as 1 and the unknown under Hd is restricted to the suspect's genotype.
Casework software such as STRmix sums over all genotypes that explain the mixture,
weighted by peak heights, so its LR will differ.
