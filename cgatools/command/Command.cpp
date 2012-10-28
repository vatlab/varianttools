// Copyright 2010 Complete Genomics, Inc.
//
// Licensed under the Apache License, Version 2.0 (the "License"); you
// may not use this file except in compliance with the License. You
// may obtain a copy of the License at
//
// http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or
// implied. See the License for the specific language governing
// permissions and limitations under the License.

#include "cgatools/core.hpp"
#include "cgatools/command/Command.hpp"
#include "cgatools/util/Exception.hpp"
#include "cgatools/util/Streams.hpp"
#include "cgatools/util/DelimitedFile.hpp"
#include <iostream>
#include <sstream>
#include <boost/algorithm/string.hpp>
#include <boost/foreach.hpp>
#include <boost/program_options/options_description.hpp>

namespace cgatools { namespace command {

    using namespace std;

    using util::Exception;
    using util::InputStream;
    using util::OutputStream;

    using boost::shared_ptr;

    namespace clstyle = po::command_line_style;
    namespace ba = boost::algorithm;

    namespace {

        const char* CGA_BETA_COMMANDS[] =
        {
            "listvariants",
            "testvariants",
            "evidence2sam",
            "evidence2cache",
            "mergedmap2sam",
            "join",
            "junctiondiff",
            "junctioncmp",
            "junctions2events",
            "generatemastervar",
            "varfilter",
            "mkvcf"
        };

    } // anonymous

    Command::Command(const std::string& name,
                     const std::string& shortDescription,
                     const std::string& formatVersionSupported,
                     const std::string& longDescription)
        : name_(name),
          formatVersionSupported_(formatVersionSupported),
          shortDescription_(shortDescription),
          longDescription_(longDescription),
          beta_(false),
          argc_(0),
          argv_(0),
#if BOOST_VERSION >= 104200
          options_(getHelpLineLength(), getHelpLineLength()-8),
#else
          options_(getHelpLineLength()),
#endif
          positionalOptions_()
    {
        options_.add_options()
            ("help,h", "Print this help message.");
        hiddenOptions_.add_options()
            ("pipeline-version-string",
             po::value<std::string>(&util::DelimitedFileMetadata::PIPELINE_VERSION),
             "pipeline software version (SOFTWARE_VERSION in metadata) override");

        BOOST_FOREACH(const char* bc, CGA_BETA_COMMANDS)
        {
            if (bc == name)
            {
                beta_ = true;
                break;
            }
        }

        if (beta_)
            options_.add_options()
                ("beta", "This is a beta command. To run this command, you must pass the --beta flag.");
    }

    std::istream& Command::openStdin (const std::string& path)
    {
        if ("STDIN" == path)
        {
            myIn_.reset(static_cast<std::istream*>(0));
            return cin;
        }
        myIn_ = InputStream::openCompressedInputStreamByExtension(path);
        return *myIn_;
    }

    std::ostream& Command::openStdout(const std::string& path)
    {
        if ("STDOUT" == path)
        {
            myOut_.reset(static_cast<std::ostream*>(0));
            return cout;
        }
        myOut_ = OutputStream::openCompressedOutputStreamByExtension(path);
        return *myOut_;
    }

    int Command::operator()(int argc, char** argv)
    {
        argc_ = argc;
        argv_ = argv;

        // Combine public options with the hidden options
        po::options_description allOptions;
        allOptions.add(options_).add(hiddenOptions_);

        po::variables_map vm;
        po::store(po::command_line_parser(argc-1, argv+1).
                  // allow_guessing is buggy -- disable it
                  style(clstyle::allow_short |
                        clstyle::short_allow_adjacent |
                        clstyle::short_allow_next |
                        clstyle::allow_long |
                        clstyle::long_allow_adjacent |
                        clstyle::long_allow_next |
                        clstyle::allow_sticky |
                        /* clstyle::allow_guessing | */
                        clstyle::allow_dash_for_short).
                  options(allOptions).
                  positional(positionalOptions_).
                  run(), vm);

        if (vm.count("help") != 0)
        {
            help(cout, false);
            return 0;
        }
        if (0 == vm.count("beta") && beta_)
            throw util::Exception("--beta flag required");

        po::notify(vm);

        return run(vm);
    }

    void Command::help(std::ostream& out, bool html)
    {
        if (html)
            out << "<a name=\"" << name_ << "\">";
        out << "COMMAND NAME";
        if (html)
            out << "</a>";
        out << "\n    " << name_ << " - " << shortDescription_ << "\n";
        if (longDescription_ != "")
            out << "\nDESCRIPTION\n" << formatDescription(longDescription_, 4, getHelpLineLength(), html);
        out << "\nOPTIONS\n";
        BOOST_FOREACH(const shared_ptr<po::option_description>& option, options_.options())
        {
            out << "  " << option->format_name() << " " << option->format_parameter() << "\n";
            out << formatDescription(option->description(), 8, getHelpLineLength(), html) << "\n";
        }
        if (formatVersionSupported_ != "")
        {
            out << "\nSUPPORTED FORMAT_VERSION\n";
            out << formatDescription(formatVersionSupported_, 4, getHelpLineLength(), html);
        }
        out << endl;
    }

    std::string Command::getCommandLine() const
    {
        string result;
        for(int ii=0; ii<argc_; ii++)
        {
            if (0 != ii)
                result += " ";
            result.push_back('\"');
            result += escapeArgument(argv_[ii]);
            result.push_back('\"');
        }
        return result;
    }

    std::string Command::escapeArgument(const std::string& arg) const
    {
        string result;
        for(size_t ii=0; ii<arg.size(); ii++)
        {
            if ('\"' == arg[ii] || '$' == arg[ii] || '\\' == arg[ii])
                result.push_back('\\');
            result.push_back(arg[ii]);
        }
        return result;
    }

    void Command::requireParam(const po::variables_map& vm, const std::string& param)
    {
        if (0 == vm.count(param))
            throw util::Exception("missing required param: "+param);
    }

    class DocToken
    {
    public:
        DocToken(const std::string& html, bool tab)
            : html_(html),
              tab_(tab)
        {
            // Only tabs can be 0-length. And tabs must be 0-length.
            CGA_ASSERT( (0 == html_.size()) == tab_ );
        }

        std::string getText() const;

        const std::string& getHtml() const
        {
            return html_;
        }

        bool isTab() const
        {
            return tab_;
        }

        bool isSpace() const
        {
            if (0 == html_.size())
                return false;
            for(size_t ii=0; ii<html_.size(); ii++)
            {
                if (' ' != html_[ii])
                    return false;
            }
            return true;
        }

        size_t textLength() const
        {
            return getText().size();
        }

        static void tokenize(std::vector<DocToken>& tokens, const std::string& line);

    private:
        static size_t addToken(std::vector<DocToken>& tokens, const std::string& line, size_t pos);
        static size_t readEntity(const std::string& line, size_t pos);

        std::string html_;
        bool tab_;
    };

    std::string DocToken::getText() const
    {
        string result;
        bool inTag = false;

        for(size_t pos=0; pos<html_.size(); pos++)
        {
            bool finishedRead = false;
            switch(html_[pos])
            {
            case '<':
                CGA_ASSERT(!inTag);
                inTag = true;
                break;
            case '>':
                CGA_ASSERT(inTag);
                inTag = false;
                finishedRead = true;
                break;
            case '&':
            {
                size_t posEnd = readEntity(html_, pos);
                string ent = ba::to_lower_copy(html_.substr(pos+1, posEnd-pos-1));
                if (!inTag)
                {
                    if ("amp" == ent)
                        result.push_back('&');
                    else if ("gt" == ent)
                        result.push_back('>');
                    else if ("lt" == ent)
                        result.push_back('<');
                    else if ("quot" == ent)
                        result.push_back('"');
                    else if ("apos" == ent)
                        result.push_back('\'');
                    else
                        throw Exception("unrecognized entity: " + html_);
                }
                pos = posEnd;
                finishedRead = true;
                break;
            }
            default:
                break;
            }

            if (finishedRead || inTag)
                continue;

            result.push_back(html_[pos]);
        }

        return result;
    }

    void DocToken::tokenize(std::vector<DocToken>& tokens, const std::string& line)
    {
        tokens.clear();
        for(size_t pos=0; pos<line.size(); pos=addToken(tokens, line, pos))
            ;
    }

    size_t DocToken::addToken(std::vector<DocToken>& tokens, const std::string& line, size_t pos)
    {
        if ('\t' == line[pos])
        {
            tokens.push_back(DocToken("", true));
            return pos+1;
        }
        if (' ' == line[pos])
        {
            size_t posEnd;
            for(posEnd=pos; posEnd<line.size(); posEnd++)
            {
                if (' ' != line[posEnd])
                    break;
            }
            tokens.push_back(DocToken(line.substr(pos, posEnd-pos), false));
            return posEnd;
        }

        // Read in a regular token.
        bool inTag = false;
        for(size_t posEnd=pos; posEnd<line.size(); posEnd++)
        {
            switch(line[posEnd])
            {
            case ' ':
                if (inTag)
                    break;
                tokens.push_back(DocToken(line.substr(pos, posEnd-pos), false));
                return posEnd;
            case '\t':
                if (inTag)
                    throw Exception("can't tokenize doc line: tab in tag: "+ line);
                tokens.push_back(DocToken(line.substr(pos, posEnd-pos), false));
                return posEnd;
            case '<':
                if (inTag)
                    throw Exception("can't tokenize doc line: < in tag: "+ line);
                inTag = true;
                break;
            case '>':
                if (!inTag)
                    throw Exception("can't tokenize doc line: > not in tag: "+ line);
                inTag = false;
                break;
            case '&':
                posEnd = readEntity(line, posEnd) + 1;
                break;
            default:
                break;
            }
        }
        tokens.push_back(DocToken(line.substr(pos, line.size()-pos), false));
        return line.size();
    }

    size_t DocToken::readEntity(const std::string& line, size_t pos)
    {
        CGA_ASSERT(pos < line.size() && '&' == line[pos]);
        for(size_t posEnd=pos; posEnd<line.size(); posEnd++)
        {
            if (';' == line[posEnd])
                return posEnd;
        }
        throw Exception("can't tokenize doc entity: invalid entity: "+ line);
    }

    std::string Command::formatLine(
        const std::string& line, size_t indent, size_t maxLineLength, bool html)
    {
        vector<DocToken> tokens;
        DocToken::tokenize(tokens, line);

        size_t tabCount = 0;
        size_t indent2 = indent;
        size_t textEaten = 0;
        BOOST_FOREACH(const DocToken& token, tokens)
        {
            if (token.isTab())
            {
                tabCount++;
                indent2 = indent + textEaten;
            }
            textEaten += token.textLength();
        }
        if (tabCount > 1)
            throw Exception("bad doc line: multiple tabs: "+ line);

        string result(indent, ' ');
        size_t lineLength = indent;
        BOOST_FOREACH(const DocToken& token, tokens)
        {
            if (lineLength + token.textLength() > maxLineLength)
            {
                result += "\n";
                result += string(indent2, ' ');
                lineLength = indent2;
                if (token.isSpace())
                    continue;
            }

            result += html ? token.getHtml() : token.getText();
            lineLength += token.textLength();
        }
        result += "\n";
        return result;
    }

    size_t Command::getHelpLineLength()
    {
        return 79;
    }

    std::string Command::formatDescription(
        const std::string& description, size_t indent, size_t maxLineLength, bool html)
    {
        vector<string> lines;
        ba::split(lines, description, ba::is_any_of("\n"));

        string result;
        BOOST_FOREACH(const string& line, lines)
            result += formatLine(line, indent, maxLineLength, html);

        return result;
    }

} } // cgatools::command
