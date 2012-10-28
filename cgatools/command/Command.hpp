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

#ifndef CGA_TOOLS_COMMAND_COMMAND_HPP_
#define CGA_TOOLS_COMMAND_COMMAND_HPP_ 1

//! @file Command.hpp
//! File containing definition of the Command class.

#include "cgatools/core.hpp"
#include <iosfwd>
#include <string>
#include <boost/program_options.hpp>

namespace cgatools { namespace command {

    namespace po = boost::program_options;

    //! The base class for writing new commands.
    class Command
    {
    public:
        //! Runs the command and returns the error code; basically,
        //! your main() is supposed to directly delegate to this
        int operator()(int argc, char** argv);

        //! Print help to out.
        void help(std::ostream& out, bool html);

        //! Return the name of this command.
        const std::string& getName() const
        {
            return name_;
        }

        //! Returns the short description of this command.
        const std::string& getShortDescription() const
        {
            return shortDescription_;
        }

        std::string getCommandLine() const;

        //! Formats description in the same way as is done for options
        //! descriptions (newline separates paragraphs, tab declares
        //! indentation on following lines of paragraph).
        //! @param indent The number of spaces to indent each line
        //!               by. This is ignored when using boost 1.41 or
        //!               before.
        //! @param lineLength The maximum number of characters in a
        //!                   line.
        //! @param html True to output html, false to output text. The
        //!             html is assumed to be within a <pre>...</pre> tag.
        static std::string formatDescription(
            const std::string& description, size_t indent, size_t maxLineLength, bool html);

        static size_t getHelpLineLength();

    protected:
        Command(const std::string& name,
                const std::string& shortDescription,
                const std::string& formatVersionSupported,
                const std::string& longDescription);

        //! Destructor can be customized, but it's preferred to complete all
        //! necessary cleanup in run() method
        virtual ~Command()
        {
        }

        //! Returns stdin if path is "STDIN", otherwise returns the
        //! result of cgatools::util::InputStream::openCompressedInputStreamByExtension().
        std::istream& openStdin (const std::string& path);

        //! Returns stdout if path is "STDOUT", otherwise returns the
        //! result of cgatools::util::OutputStream::openCompressedOutputStreamByExtension().
        std::ostream& openStdout(const std::string& path);

        //! Throws an exception if the given param was not specified on
        //! the command line.
        void requireParam(const po::variables_map& vm, const std::string& param);

        //! The actual code of the command goes here. Throw exceptions to
        //! communicate errors.
        virtual int run(po::variables_map& vm) = 0;

    private:
        std::string escapeArgument(const std::string& arg) const;
        static std::string formatLine(
            const std::string& line, size_t indent, size_t maxLineLength, bool html);

        std::string name_;
        std::string formatVersionSupported_;
        std::string shortDescription_;
        std::string longDescription_;
        bool beta_;
        boost::shared_ptr<std::istream> myIn_;
        boost::shared_ptr<std::ostream> myOut_;

        int argc_;
        char** argv_;

    protected:
        po::options_description options_;
        po::options_description hiddenOptions_;
        po::positional_options_description positionalOptions_;
    };

} } // cgatools::command

#endif // CGA_TOOLS_COMMAND_COMMAND_HPP_
