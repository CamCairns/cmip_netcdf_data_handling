## Counting front-page NYT bylines with Bash and jq

A quickie script in Bash, with the [jq JSON-command-line-parser](http://stedolan.github.io/jq/), to access the [New York Times Article Search API v2](http://developer.nytimes.com/docs/read/article_search_api_v2) and count up the bylines. You can count whatever you want obviously, but a student was interested in replicating the gender calculation found at [Who Writes for the New York Times?](http://www.whowritesfor.com/).


### The NYT Article Search API

First thing you have to do is signup and register as a developer. API Keys are assigned by API, so make sure you specify the [Article Search API](http://developer.nytimes.com/docs/read/article_search_api_v2).

Even before you register, you can use the NYT's handy API Console to interactively test your queries: http://developer.nytimes.com/io-docs

The [Article Search API](http://developer.nytimes.com/docs/read/article_search_api_v2) is pretty flexible; you can call it with no parameters except for your `api-key` and it will return (presumably) a list of articles, in reverse chronological order, starting from Sept. 18, 1851. However, it only returns 10 articles per request. And it won't let you paginate beyond a `page` parameter of `100` (i.e. you can't go to `page` 100000 to retrieve the 1,000,000th oldest Times article). To put it another way, you can only paginate through a maximum of 10,000 results, so you'll have to facet your search.

In the [shell script attached to this gist](https://gist.github.com/dannguyen/07e91763f1f5fd410c84#file-nyt-articles-by-source-collector-sh), I set a few parameters to limit the number of possible results (the API refers to them as "`hits`"):

- `begin_date` - specifies the oldest day to include in the search, and takes in a date formatted as `YYYYMMDD`, e.g. `20150316`
- `end_date` - specifies the most recent day to include in the search. If set equal to `begin_date`, you can effectively limit your search to that day
- `fq` - This parameter lets you filter results with [Lucene syntax](http://lucene.apache.org/core/2_9_4/queryparsersyntax.html). I don't know what that is so I just Googled around and [copied from this iPython notebook](http://nbviewer.ipython.org/github/Jay-Oh-eN/pydatasv2014/blob/master/example.ipynb). The API returns objects with a `source` field. I just want articles from the _New York Times_ (as opposed to wire news from _Reuters_ and _Associated Press_), so I use this key-value parameter in the API call: `fq=source.contains:("New York")`
- `page` - Increment this parameter (starting at `0` for results `0` to `9`) to paginate the results.

For any given day, there seems to be around 200 to 300 NYT-bylined articles. Here's the [first page (`page=0`) of results for Mar. 15, 2015](https://gist.github.com/dannguyen/07e91763f1f5fd410c84#file-20150316-page-0-json)

In my [sample shell script](https://gist.github.com/dannguyen/07e91763f1f5fd410c84#file-nyt-articles-by-source-collector-sh), each day is downloaded into a subdirectory such as `./data-hold/20150316`.

To get the bylines &ndash; last name, first name, and listed position (i.e. `rank`) &ndash; for a given day, after the data's been downloaded, you can use __jq__. 

The following __jq__ query uses the `select` function to filter the list to __front page__ bylines:


~~~sh
cat data-hold/20150316/*.json | \
  jq -r '.response .docs[] | 
  select(.print_page == "1") .byline .person[] | 
  [.lastname, .firstname, .rank] | @csv'
~~~

&ndash; which results in this

```
"ELIGON","John",1
"BRANTLEY","Ben",1
"SISARIO","Ben",1
"HOROWITZ","Jason",1
"RICHTEL","Matt",1
"TRACY","Marc",1
"GOLDSTEIN","Matthew",1
"PERLROTH","Nicole",2
"WEISMAN","Jonathan",1
"CORKERY","Michael",1
"SILVER-GREENBERG","Jessica",2
"HADID","Diaa",1
"DOUGHERTY","Conor",1
"HARDY","Quentin",2
"CARAMANICA","Jon",1
"GENZLINGER","Neil",1
"TABUCHI","Hiroko",1
"TRACY","Marc",1
"GRIMES","William",1
"BAGLI","Charles",1
"YEE","Vivian",2
"SIEGAL","Nina",1
```


To include the actual headlines per byline, it's kind of a pain in the ass to do it just from __jq__, so might as well use Python: 

~~~py
import json
import glob
print('|', '|'.join(('lastname', 'firstname', 'rank', 'desk', 'headline')), '|')
print('|--|--|--|--|--|')
for filename in glob.glob('./data-hold/20150316/*.json'):
    data = json.loads(open(filename).read())
    for article in data['response']['docs']:
        if article['print_page'] == "1":
            headline = article['headline']['main']
            desk = article['news_desk']        
            for p in article['byline']['person']:
                print('|', '|'.join((p['lastname'], p['firstname'], str(p['rank']), desk, headline)), '|')
~~~


The results, as Markdown-formatted tables:

|     lastname     | firstname | rank |   desk   |                                       headline                                       |
|------------------|-----------|------|----------|--------------------------------------------------------------------------------------|
| ELIGON           | John      |    1 | National | Crackdown in a Detroit Stripped of Metal Parts                                       |
| BRANTLEY         | Ben       |    1 | Culture  | Review: ‘On the Twentieth Century,’ With Kristin Chenoweth, Opens on Broadway        |
| SISARIO          | Ben       |    1 | Business | ‘Blurred Lines’ Lawyer Rocks Music Industry Again                                    |
| HOROWITZ         | Jason     |    1 | National | Evangelicals Aim to Mobilize an Army for Republicans in 2016                         |
| RICHTEL          | Matt      |    1 | Business | A Police Gadget Tracks Phones? Shhh! It&#8217;s Secret                               |
| TRACY            | Marc      |    1 | Sports   | Kentucky, No. 1 in Height, Too, Relishes the View From on High                       |
| GOLDSTEIN        | Matthew   |    1 | Business | Authorities Closing In on Hackers Who Stole Data From JPMorgan Chase                 |
| PERLROTH         | Nicole    |    2 | Business | Authorities Closing In on Hackers Who Stole Data From JPMorgan Chase                 |
| WEISMAN          | Jonathan  |    1 | National | Chasm Grows Within G.O.P. Over Spending                                              |
| CORKERY          | Michael   |    1 | Business | Many Buyers for Subprime Auto Loan Bundle                                            |
| SILVER-GREENBERG | Jessica   |    2 | Business | Many Buyers for Subprime Auto Loan Bundle                                            |
| HADID            | Diaa      |    1 | Foreign  | Arab Alliance Rises as Force in Israeli Elections                                    |
| DOUGHERTY        | Conor     |    1 | Business | Managers Turn to Computer Games, Aiming for More Efficient Employees                 |
| HARDY            | Quentin   |    2 | Business | Managers Turn to Computer Games, Aiming for More Efficient Employees                 |
| CARAMANICA       | Jon       |    1 | Culture  | Without Joan Rivers, ‘Fashion Police’ Is Falling Apart                               |
| GENZLINGER       | Neil      |    1 | Culture  | Review: ‘iZombie,’ the Undead as a Force for Good                                    |
| TABUCHI          | Hiroko    |    1 | Business | Etsy’s Success Gives Rise to Problems of Credibility and Scale                       |
| TRACY            | Marc      |    1 | Sports   | N.C.A.A. Tournament 2015: Villanova, Duke and Wisconsin Join Kentucky as No. 1 Seeds |
| GRIMES           | William   |    1 | Culture  | Asia Week Is Highlighted by Robert Hatfield Ellsworth Sale                           |
| BAGLI            | Charles   |    1 | Metro    | Robert Durst of HBO’s ‘The Jinx’ Says He ‘Killed Them All’                           |
| YEE              | Vivian    |    2 | Metro    | Robert Durst of HBO’s ‘The Jinx’ Says He ‘Killed Them All’                           |
| SIEGAL           | Nina      |    1 | Culture  | For Michaela DePrince, a Dream Comes True at the Dutch National Ballet               |