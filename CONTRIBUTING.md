# Contributing to KrakenFOAM

Thank you for your interest in KrakenFOAM. Bug reports, feature requests and
pull requests are welcome.

## Reporting issues

Please open an issue on GitHub and include:

- the OpenFOAM version (`foamVersion`) and operating system,
- the Cerberus/Kraken version used,
- a description of the problem and, if possible, a minimal case reproducing
  it together with the relevant log files.

## Pull requests

1. Fork the repository and create a branch for your change.
2. Make sure the code builds without warnings against OpenFOAM 13:

   ```sh
   source /opt/openfoam13/etc/bashrc
   ./Allwclean && ./Allwmake -j
   ```

3. Check that the tutorials affected by the change still run.
4. Keep commits focused and describe *what* and *why* in the commit message.
5. Open a pull request against the `main` branch.

## Coding style

KrakenFOAM follows the
[OpenFOAM coding style](https://openfoam.org/dev/coding-style-guide/) used in
OpenFOAM 13. In short:

- Maximum line length of 80 characters, 4-space indentation, no tabs, no
  trailing whitespace and Unix (LF) line endings.
- Opening and closing braces on their own lines; function arguments that do
  not fit on one line are split one per line inside `(` `)` on separate
  lines.
- No spaces around `*` and `/` in expressions, `>>` without a space in nested
  templates, `forAll(list, i)`, `Info<<`, `nullptr`.
- Document classes in the header banner (`Description`, `Usage`,
  `SourceFiles`) and members with `//-` comments.
- Section separators such as
  `// * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //`
  and the closing
  `// ************************************************************************* //`
  line, as in the OpenFOAM sources.
- Case dictionaries use the OpenFOAM 13 file banner.

The OpenFOAM git hooks (`$WM_PROJECT_DIR/bin/tools/pre-commit-hook`) check
most of these rules and can be used as a reference.

## License

By contributing you agree that your contributions are licensed under the
GNU General Public License v3.0 or later, see [LICENSE](LICENSE). New source
files should carry the same license header as the existing files.
