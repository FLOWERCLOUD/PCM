#include "main_window.h"
#include <QApplication>
#include <iostream>
#include <QLocale>
#include <QTranslator>
#include <QTextCodec>
#include <QDebug>
#include "console.h"

class MyApplication : public QApplication
{
public:
	MyApplication(int& argc, char ** argv) : QApplication(argc, argv) { }
	virtual ~MyApplication() { }

	// reimplemented from QApplication so we can throw exceptions in slots
	virtual bool notify(QObject * receiver, QEvent * event) {
		try {
			return QApplication::notify(receiver, event);
		} catch(std::exception& e) {
			qCritical() << "Exception thrown:" << e.what();
		}
		return false;
	}
};

int main(int argc, char *argv[])
{
	ConsoleLogger::Instance();
	MyApplication a(argc, argv);
	main_window w;
	w.show();
	return a.exec();
}
