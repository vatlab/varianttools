from pyparsing_py2 import *

condition_items = Forward().setName("conditions in where statement")
result_items = Forward().setName("items in select statement")

LPAR,RPAR,COMMA = "(", ")", ","
# define all keywords
(UNION, ALL, AND, INTERSECT, EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, 
 CROSS, LEFT, OUTER, JOIN, AS, INDEXED, NOT, SELECT, DISTINCT, FROM, WHERE, GROUP, BY,
 HAVING, ORDER, BY, LIMIT, OFFSET, CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, END,
 CASE, WHEN, THEN, EXISTS, COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, CURRENT_TIME,
 CURRENT_DATE, CURRENT_TIMESTAMP) = map(CaselessKeyword, """UNION, ALL, AND, INTERSECT, 
 EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, CROSS, LEFT, OUTER, JOIN, AS,
 INDEXED, NOT, SELECT, DISTINCT, FROM, WHERE, GROUP, BY, HAVING, ORDER, BY, LIMIT,
 OFFSET, CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, END, CASE, WHEN, THEN, EXISTS,
 COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, CURRENT_TIME, CURRENT_DATE,
 CURRENT_TIMESTAMP""".replace(",","").split())

keyword = MatchFirst((UNION, ALL, INTERSECT, EXCEPT, COLLATE, ASC, DESC, ON, USING, NATURAL, INNER, 
 CROSS, LEFT, OUTER, JOIN, AS, INDEXED, NOT, SELECT, DISTINCT, FROM, WHERE, GROUP, BY,
 HAVING, ORDER, BY, LIMIT, OFFSET, CAST, ISNULL, NOTNULL, NULL, IS, BETWEEN, ELSE, END, CASE, WHEN,
 THEN, EXISTS, COLLATE, IN, LIKE, GLOB, REGEXP, MATCH, ESCAPE, CURRENT_TIME, CURRENT_DATE, 
 CURRENT_TIMESTAMP))

identifier = ~keyword + Word(alphas, alphanums+"_")
table_name = identifier.copy().setName('table_name')

# * or table.*
all_fields = Group(Optional(table_name + ".") + Literal('*')).setName('all fields')
# sift_score or dbNSFP.sift_score
field = Group(Optional(table_name + ".") + identifier).setName('field')
# count(1)
integer = Regex(r"[+-]?\d+")
numeric_literal = Regex(r"\d+(\.\d*)?([eE][+-]?\d+)?")
string_literal = QuotedString("'") | QuotedString('"')
literal_value = (numeric_literal | string_literal | NULL | CURRENT_TIME | CURRENT_DATE | CURRENT_TIMESTAMP )
expr_term = Group(identifier + LPAR + identifier + RPAR | \
    identifier + LPAR + Optional(delimitedList(identifier)) + RPAR | \
    literal_value)
# pos - 1, pos is not NULL, func(pos)

UNARY,BINARY,TERNARY=1,2,3
expr = Group(operatorPrecedence(expr_term,
    [
    (oneOf('- + ~') | NOT, UNARY, opAssoc.LEFT),
    ('||', BINARY, opAssoc.LEFT),
    (oneOf('* / %'), BINARY, opAssoc.LEFT),
    (oneOf('+ -'), BINARY, opAssoc.LEFT),
    (oneOf('<< >> & |'), BINARY, opAssoc.LEFT),
    (oneOf('< <= > >='), BINARY, opAssoc.LEFT),
    (oneOf('= == != <>') | IS | IN | LIKE | GLOB | MATCH | REGEXP, BINARY, opAssoc.LEFT),
    ('||', BINARY, opAssoc.LEFT),
    ((BETWEEN,AND), TERNARY, opAssoc.LEFT),
    ]))

result_column = all_fields | field | expr_term


result_items << delimitedList(result_column)

fields = []

import sys
if __name__ == '__main__':
    print sys.argv[1]
    try:
        #print result_items.parseString(sys.argv[1]).dump()
        print result_items.parseString(sys.argv[1])
        print fields
    except ParseException, pe:
        print pe.msg
    print
