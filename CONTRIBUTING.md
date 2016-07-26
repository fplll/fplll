# Contributing to fplll #

fplll welcomes contributions. We encourage users to fix bugs, improve the documentation, write tests and to enhance the code. We encourage researchers to contribute implementations of their algorithms to fplll. In the following we are trying to give some guidance on how to contribute effectively.

## Communication ##

Communication in the fplll project happens mainly on [GitHub](https://github.com/fplll/fplll/issues) and on our [mailing list](https://groups.google.com/forum/#!forum/fplll-devel).

## Setup ##

### Fork on GitHub ##

Login/signup on GitHub and fork fplll from [https://github.com/fplll/fplll](https://github.com/fplll/fplll).

### Clone your fork locally

Now clone your git repo where my-github-name is your account name on GitHub:

    $ git clone git@github.com:my-github-name/fplll.git

Then run

    $ ./autogen.sh
    $ ./configure
    $ make
    $ make check

as usual.

## Reporting Bugs ##

Bug should be filed at [https://github.com/fplll/fplll/issues](https://github.com/fplll/fplll/issues). Alternatively, feel free to contact  [https://groups.google.com/forum/#!forum/fplll-devel](https://groups.google.com/forum/#!forum/fplll-devel). The former method is definitely preferred though.

## Setting up topic branches and generating pull requests

While it's handy to provide useful code snippets in an issue, it is better for to submit pull requests.

In Git it is best to isolate each topic or feature into a "topic branch". While individual commits allow you control over how small individual changes are made to the code, branches are a great way to group a set of commits all related to one feature together, or to isolate different efforts when you might be working on multiple topics at the same time.

While it takes some experience to get the right feel about how to break up commits, a topic branch should be limited in scope to a single ``issue`` as submitted to an issue tracker.

Also since GitHub pegs and syncs a pull request to a specific branch, it is the **only** way that you can submit more than one fix at a time. If you submit a pull from your ``master`` branch, you can't make any more commits to your ``master`` branch without those getting added to the pull.

To create a topic branch, its easiest to use the convenient ``-b`` argument to ``git checkout``

    $ git checkout -b fix-broken-thing
    Switched to a new branch 'fix-broken-thing'

You should use a verbose enough name for your branch so it is clear what it is about. Now you can commit your changes and regularly merge in the upstream develop as described below.

When you are ready to generate a pull request, either for preliminary review, or for consideration of merging into the project you must first push your local topic branch back up to GitHub

    git push origin fix-broken-thing

Now when you go to your fork on GitHub, you will see this branch listed under the "Source" tab where it says "Switch Branches". Go ahead and select your topic branch from this list, and then click the "Pull request" button.

Here you can add a comment about your branch. If this in response to a submitted issue, it is good to put a link to that issue in this initial comment. The maintainers will be notified of your pull request and it will be reviewed (see below for best practices). Note that you can continue to add commits to your topic branch (and push them up to GitHub) either if you see something that needs changing, or in response to a reviewer's comments. If a reviewer asks for changes, you do not need to close the pull and reissue it after making changes. Just make the changes locally, push them to GitHub, then add a comment to the discussion section of the pull request.

## Pull upstream changes into your fork regularly

It is important that you pull upstream changes from master into your fork on a regular basis. Nothing is worse than putting in a days of hard work into a pull request only to have it rejected because it has diverged too far from master. 

To pull in upstream changes::

    $ git remote add upstream https://github.com/fplll/fplll.git
    $ git fetch upstream master

Check the log to be sure that you actually want the changes, before merging::

    $ git log upstream/master

Then merge the changes that you fetched::

    $ git merge upstream/master

For more info, see http://help.github.com/fork-a-repo/

## How to get your pull request accepted

We want your submission. But we also want to provide a stable experience for our users. In the following, we give some guidelines for contributing.

### Run the tests!

Before you submit a pull request, please run tests:

    $ make check

These checks are also run on [Travis-CI](https://travis-ci.org/fplll/fplll) automatically for every pull request. Nothing failing tests will be accepted.

### If you add code please add tests

Code that isn’t tested is broken.

Also, keep your tests simple. Complex tests end up requiring their own tests. We would rather see duplicated assertions across test methods then cunning utility methods that magically determine which assertions are needed at a particular stage. Remember: Explicit is better than implicit.

Finally, the nature of fplll means that sometimes it is hard to properly test the behaviour of a change quickly. Running BKZ for several minutes takes way too long for a test. In this case, we should at least test that a particular piece of code compiles and runs.

### Keep your pull requests limited to a single issue

Pull requests should be as small/atomic as possible.

### Coding Conventions

fplll is written in [C++11](https://en.wikipedia.org/wiki/C%2B%2B11) and we try to make use of its modern features to make the library readable.

Please keep your code as clean and straightforward as possible. Code is written for the consumption by compilers and for the consumption by human beings. By making code clear and easy to understand, others can build on it and fix issues should they arise.

Our naming convention is close to Python's [naming convention](https://www.python.org/dev/peps/pep-0008/). Classes are in ``CamelCase``. Functions, methods, parameters and local variables in ``lower_case`` . Curly braces go on the next line and we [prefer explicit curly braces](https://nakedsecurity.sophos.com/2014/02/24/anatomy-of-a-goto-fail-apples-ssl-bug-explained-plus-an-unofficial-patch/), e.g.

```c++
if (foo)
{
  do_something_good();
}
```
    
instead of:

```c++
if (foo)
  do_something_bad();
```

The following [clang-format](http://clang.llvm.org/docs/ClangFormat.html) config might help to format your code.

```yaml
---
Language:        Cpp
# BasedOnStyle:  LLVM
AccessModifierOffset: -2
AlignAfterOpenBracket: true
AlignConsecutiveAssignments: true
AlignEscapedNewlinesLeft: false
AlignOperands:   true
AlignTrailingComments: true
AllowAllParametersOfDeclarationOnNextLine: true
AllowShortBlocksOnASingleLine: false
AllowShortCaseLabelsOnASingleLine: false
AllowShortIfStatementsOnASingleLine: false
AllowShortLoopsOnASingleLine: false
AllowShortFunctionsOnASingleLine: All
AlwaysBreakAfterDefinitionReturnType: false
AlwaysBreakTemplateDeclarations: false
AlwaysBreakBeforeMultilineStrings: false
BreakBeforeBinaryOperators: None
BreakBeforeTernaryOperators: true
BreakConstructorInitializersBeforeComma: false
BinPackParameters: true
BinPackArguments: true
ColumnLimit:     100
ConstructorInitializerAllOnOneLineOrOnePerLine: false
ConstructorInitializerIndentWidth: 4
DerivePointerAlignment: false
ExperimentalAutoDetectBinPacking: false
IndentCaseLabels: false
IndentWrappedFunctionNames: false
IndentFunctionDeclarationAfterType: false
MaxEmptyLinesToKeep: 1
KeepEmptyLinesAtTheStartOfBlocks: true
NamespaceIndentation: None
ObjCBlockIndentWidth: 2
ObjCSpaceAfterProperty: false
ObjCSpaceBeforeProtocolList: true
PenaltyBreakBeforeFirstCallParameter: 19
PenaltyBreakComment: 300
PenaltyBreakString: 1000
PenaltyBreakFirstLessLess: 120
PenaltyExcessCharacter: 1000000
PenaltyReturnTypeOnItsOwnLine: 60
PointerAlignment: Right
SpacesBeforeTrailingComments: 2
Standard:        Cpp11
IndentWidth:     2
TabWidth:        8
UseTab:          Never
BreakBeforeBraces: Allman
SpacesInParentheses: false
SpacesInSquareBrackets: false
SpacesInAngles:  false
SpaceInEmptyParentheses: false
SpacesInCStyleCastParentheses: false
SpaceAfterCStyleCast: false
SpacesInContainerLiterals: true
SpaceBeforeAssignmentOperators: true
ContinuationIndentWidth: 4
CommentPragmas:  '^ IWYU pragma:'
ForEachMacros:   [ foreach, Q_FOREACH, BOOST_FOREACH ]
SpaceBeforeParens: ControlStatements
DisableFormat:   false
...
```


Furthermore, the pixel shortage is over. We want to see:

- `package` instead of `pkg`
- `grid` instead of `g`
- `my_function_that_does_things` instead of `mftdt`

### Faster compilation

By default, libtool builds everything twice, one for the static and one for the dynamic library, cf. https://stackoverflow.com/questions/572760/libtool-slowness-double-building. If you want to avoid this double compiling time you can run ./configure --disable-static which disables building the static library.

### Attribution

Please do not forget to add yourself as a contributor in [README.md](README.md) if you make a non-trivial contribution. Furthermore, you may want to claim copyright in the copyright headers of each file.

## Documentation

fplll uses [doxygen](http://www.stack.nl/~dimitri/doxygen/) with a [bootstrap theme](https://github.com/Velron/doxygen-bootstrapped) to generate API documentation. To produce API documentation run

    $ doxygen Doxyfile

Our documentation is served at [https://fplll.github.io/fplll/](https://fplll.github.io/fplll/) using [GitHub pages](https://pages.github.com). To update the documentation, check out the ``gh-pages`` branch and update the html files in there. Doxygen writes its outputs to ``doc/html``, you can arrange it that this directory holds the ``gh-pages`` branch of the fplll repository:

    $ cd doc
    $ git clone -b gh-pages git@github.com::<my-github-name>/fplll.git html
    $ cd ..
    
Now, whenever you run ``doxygen`` it will write its outputs to a directory which holds the right branch. If you push it to your remote, you can then check it at [http://my-github-name.github.io/fplll](http://my-github-name.github.io/fplll).
