import QtQuick 2.15
import QtQuick.Controls 2.15

ToolButton {
    id: toolButton
    text: qsTr("Tool Button")
    property string pageSource: ;
    property int pageIndex: 1
    onClicked: {
        //console.log(Qt.application.engine.style)
        loader.source = pageSource
        swipeView.currentIndex = pageIndex
    }
}
