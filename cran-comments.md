This is a re-submission.

## Test environments
* Ubuntu 16.04 (travis): r-devel (2020-03-13 r77948), r-release (3.6.2), r-oldrelease (3.5.3)

* winbuilder: r-devel (2020-03-11 r77925), r-release (3.6.3)

## R CMD check results
> There were no ERRORs or WARNINGs.

> There was 1 NOTE:

> Maintainer: ‘Alexander Rix <alexrix@umich.edu>’

> New submission

## CRAN Reviewer Comments

> Thanks, your examples are wrapped in \donttest{}, hence nothing gets
> tested. Please unwrap the examples if that is feasible and if they can
> be executed in < 5 sec for each Rd file or create additionally small toy
> examples. Something like ... Alternatively, you can write some tests
> (e.g. by using testthat). The overall execution time should not exceed 10 minutes.

I added the examples as `testhat` tests. These tests take ~250s in total for
winbuilder.

> Missing Rd-tags:
>      coef.cv.galasso.Rd: \value
>      coef.cv.saenet.Rd: \value
>      coef.galasso.Rd: \value
>      coef.saenet.Rd: \value
>
> Please add the tag and explain in detail the returned objects.

This has been fixed
